from typing import TextIO, Dict, Tuple
from collections import defaultdict

from .loader import open_file, load_barcode_whitelist, load_barcode_remapping
from .matcher import find_closest_barcode, remap_barcode

# Translation table for reverse complement — faster than dict lookup per base
_RC_TABLE = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
_BASES = 'ACGT'


def _match_barcode(bc: str, whitelist_set: set, max_mismatches: int, quality: str):
    """
    Match bc against whitelist using set membership — O(1) exact, O(48) for
    1-mismatch, O(1080) for 2-mismatches — instead of O(737K) linear scan.

    Returns (matched_barcode, num_mismatches) or (None, None).
    """
    # Exact match
    if bc in whitelist_set:
        return bc, 0

    if max_mismatches < 1:
        return None, None

    # 1-mismatch: generate all 16×3 = 48 variants and probe the set
    hits_1 = []
    for i in range(len(bc)):
        orig = bc[i]
        for base in _BASES:
            if base != orig:
                variant = bc[:i] + base + bc[i+1:]
                if variant in whitelist_set:
                    hits_1.append((variant, i))

    if hits_1:
        if len(hits_1) == 1 or not quality or len(quality) != len(bc):
            return hits_1[0][0], 1
        # Tiebreak: prefer lower quality score at mismatch position (less confident base)
        return min(hits_1, key=lambda h: ord(quality[h[1]]))[0], 1

    if max_mismatches < 2:
        return None, None

    # 2-mismatch: C(16,2)×9 = 1080 variants
    hits_2 = []
    for i in range(len(bc)):
        for j in range(i + 1, len(bc)):
            for b1 in _BASES:
                if b1 == bc[i]:
                    continue
                for b2 in _BASES:
                    if b2 == bc[j]:
                        continue
                    variant = bc[:i] + b1 + bc[i+1:j] + b2 + bc[j+1:]
                    if variant in whitelist_set:
                        hits_2.append((variant, i, j))

    if not hits_2:
        return None, None
    if len(hits_2) == 1 or not quality or len(quality) != len(bc):
        return hits_2[0][0], 2
    return min(hits_2, key=lambda h: ord(quality[h[1]]) + ord(quality[h[2]]))[0], 2


class BarcodeCorrector:

    def __init__(
        self,
        dna_whitelist_file: str,
        rna_whitelist_file: str,
        max_mismatches: int = 1,
        min_frac_bcs_to_find: float = 0.5,
        barcode_suffix: str = "1"
    ):
        self.whitelist_barcodes = load_barcode_whitelist(dna_whitelist_file)
        self.remapping_dict = load_barcode_remapping(dna_whitelist_file, rna_whitelist_file)
        self.max_mismatches = max_mismatches
        self.min_frac_bcs_to_find = min_frac_bcs_to_find
        self.barcode_suffix = barcode_suffix
        self.total_reads = 0
        self.reads_with_correctable_bc = 0
        self.reads_with_remapped_bc = 0
        self.mismatch_distribution = defaultdict(int)

    @staticmethod
    def reverse_complement(sequence: str) -> str:
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, 'N') for base in reversed(sequence.upper()))

    def correct_barcode(self, raw_barcode: str, quality: str = None) -> Tuple[str, bool]:
        closest_barcode = find_closest_barcode(
            raw_barcode,
            self.whitelist_barcodes,
            self.max_mismatches,
            quality=quality
        )

        if closest_barcode is None:
            return None, False

        mismatches = sum(c1 != c2 for c1, c2 in zip(raw_barcode, closest_barcode))
        self.mismatch_distribution[mismatches] += 1

        remapped_barcode = remap_barcode(closest_barcode, self.remapping_dict)
        if remapped_barcode != closest_barcode:
            self.reads_with_remapped_bc += 1

        return remapped_barcode, True

    def process_fastq(
        self,
        input_fastq: str,
        output_file: TextIO,
        num_threads: int = 1
    ) -> None:
        whitelist_barcodes = self.whitelist_barcodes
        remapping_dict = self.remapping_dict
        max_mismatches = self.max_mismatches
        barcode_suffix = self.barcode_suffix

        with open_file(input_fastq) as fh:
            while True:
                name_line = fh.readline()
                if not name_line:
                    break

                seq_line  = fh.readline().strip()
                fh.readline()                          # + separator
                qual_line = fh.readline().strip()

                read_name = name_line.strip()
                if read_name.startswith('@'):
                    read_name = read_name[1:]

                self.total_reads += 1

                seq_no_sep  = seq_line[1:]  if seq_line.startswith('N')  else seq_line
                qual_no_sep = qual_line[1:] if len(qual_line) == len(seq_line) and seq_line.startswith('N') else qual_line

                right_16_bases = seq_no_sep[-16:]
                right_16_qual  = qual_no_sep[-16:]

                # Reverse complement via fast translate table, no per-base dict lookup
                bc_for_matching = right_16_bases[::-1].translate(_RC_TABLE).upper()

                matched_bc, num_mismatches = _match_barcode(
                    bc_for_matching, whitelist_barcodes, max_mismatches, right_16_qual
                )

                cr_tag = bc_for_matching
                cy_tag = right_16_qual[::-1]

                if matched_bc is not None:
                    self.reads_with_correctable_bc += 1
                    self.mismatch_distribution[num_mismatches] += 1
                    remapped = remapping_dict.get(matched_bc, matched_bc)
                    if remapped != matched_bc:
                        self.reads_with_remapped_bc += 1
                    output_file.write(f"{read_name} CR:Z:{cr_tag}\tCY:Z:{cy_tag}\tCB:Z:{remapped}-{barcode_suffix}\n")
                else:
                    output_file.write(f"{read_name} CR:Z:{cr_tag}\tCY:Z:{cy_tag}\n")

    def get_statistics(self) -> Dict:
        frac_corrected = (
            self.reads_with_correctable_bc / self.total_reads
            if self.total_reads > 0 else 0
        )

        return {
            'total_reads': self.total_reads,
            'reads_with_correctable_bc': self.reads_with_correctable_bc,
            'fraction_correctable': frac_corrected,
            'reads_with_remapped_bc': self.reads_with_remapped_bc,
            'mismatch_distribution': dict(self.mismatch_distribution),
        }

    def write_statistics(self, stats_file: str) -> None:
        stats = self.get_statistics()
        with open(stats_file, 'w') as fh:
            fh.write("metric\tvalue\n")
            fh.write(f"total_reads\t{stats['total_reads']}\n")
            fh.write(f"reads_with_correctable_bc\t{stats['reads_with_correctable_bc']}\n")
            fh.write(f"fraction_correctable\t{stats['fraction_correctable']:.4f}\n")
            fh.write(f"reads_with_remapped_bc\t{stats['reads_with_remapped_bc']}\n")

            for mismatches in sorted(stats['mismatch_distribution'].keys()):
                count = stats['mismatch_distribution'][mismatches]
                fh.write(f"mismatches_{mismatches}\t{count}\n")

import gzip
from pathlib import Path
from typing import TextIO, Dict, Tuple
from collections import defaultdict

from .loader import open_file, load_barcode_whitelist, load_barcode_remapping
from .matcher import find_closest_barcode, remap_barcode


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
        output_file: TextIO
    ) -> None:
        with open_file(input_fastq) as fh:
            while True:
                read_name_line = fh.readline().strip()
                if not read_name_line:
                    break

                sequence_line = fh.readline().strip()
                plus_line = fh.readline().strip()
                quality_line = fh.readline().strip()

                read_name = read_name_line[1:] if read_name_line.startswith('@') else read_name_line
                raw_bc = sequence_line
                raw_bc_qual = quality_line

                self.total_reads += 1

                seq_no_sep = raw_bc[1:] if raw_bc.startswith('N') else raw_bc
                qual_no_sep = raw_bc_qual[1:] if len(raw_bc_qual) == len(raw_bc) and raw_bc.startswith('N') else raw_bc_qual
                
                right_16_bases = seq_no_sep[-16:] if len(seq_no_sep) >= 16 else seq_no_sep
                right_16_qual = qual_no_sep[-16:] if len(qual_no_sep) >= 16 else qual_no_sep
                
                bc_for_matching = self.reverse_complement(right_16_bases)
                corrected_bc, was_corrected = self.correct_barcode(bc_for_matching, quality=right_16_qual)
                
                cr_tag = bc_for_matching
                cy_tag = right_16_qual[::-1]

                if was_corrected:
                    self.reads_with_correctable_bc += 1
                    corrected_bc_with_suffix = f"{corrected_bc}-{self.barcode_suffix}"
                    output_file.write(
                        f"{read_name} CR:Z:{cr_tag}\tCY:Z:{cy_tag}\tCB:Z:{corrected_bc_with_suffix}\n"
                    )
                else:
                    output_file.write(
                        f"{read_name} CR:Z:{cr_tag}\tCY:Z:{cy_tag}\n"
                    )

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

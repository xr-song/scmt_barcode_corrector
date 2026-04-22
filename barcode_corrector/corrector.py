from typing import TextIO, Dict, Tuple, Iterator, List
import gzip
import os
import shutil

from collections import defaultdict
from multiprocessing import Pool

from .loader import open_file, load_barcode_whitelist, load_barcode_remapping
from .matcher import find_closest_barcode, remap_barcode

# Translation table for reverse complement — faster than dict lookup per base
_RC_TABLE = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
_BASES = 'ACGT'

# Per-worker process globals — set once at worker startup via Pool initializer,
# avoiding re-pickling the large whitelist set on every task.
_wl: set = None
_remap: dict = None
_max_mm: int = None
_suffix: str = None


def _worker_init(whitelist, remapping, max_mismatches, barcode_suffix):
    global _wl, _remap, _max_mm, _suffix
    _wl = whitelist
    _remap = remapping
    _max_mm = max_mismatches
    _suffix = barcode_suffix


def _match_barcode(bc: str, whitelist_set: set, max_mismatches: int, quality: str):
    """
    Match bc against whitelist using set membership — O(1) exact, O(48) for
    1-mismatch, O(1080) for 2-mismatches — instead of O(737K) linear scan.

    `quality` must be in the same coordinate space as `bc` (i.e. already
    reversed if bc is a reverse complement). Position i in quality corresponds
    to position i in bc.

    Tie handling:
      - Single hit: returned immediately.
      - Multiple hits at equal distance: prefer the candidate whose mismatch
        falls at the position with the LOWEST quality score in the read
        (low Phred = low confidence = most likely a sequencing error there).
      - If quality is absent or mismatched in length: first hit is returned
        (non-deterministic across runs due to set ordering).

    Returns (matched_barcode, num_mismatches) or (None, None).
    """
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
        # quality[i] is the Phred score of position i in bc — lower = less confident
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


def _read_chunks(fh, chunk_size: int) -> Iterator[List[Tuple[str, str, str]]]:
    """Lazy FASTQ iterator: yields fixed-size record lists."""
    chunk = []
    while True:
        name_line = fh.readline()
        if not name_line:
            if chunk:
                yield chunk
            break
        seq_line  = fh.readline().strip()
        fh.readline()  # + separator
        qual_line = fh.readline().strip()
        chunk.append((name_line.strip().lstrip('@'), seq_line, qual_line))
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []


def _split_fastq(input_fastq: str, n_parts: int, tmp_dir: str) -> List[str]:
    """
    Split input_fastq into n_parts sequential gzipped chunk files under tmp_dir.
    First pass counts total records; second pass writes equal-sized blocks so
    that concatenating chunk outputs in order reproduces the original read order.
    Returns list of paths for non-empty chunks.
    """
    # Count records (4 lines each) with a fast line count
    n_records = 0
    with open_file(input_fastq) as fh:
        for _ in fh:
            n_records += 1
    n_records //= 4

    records_per_chunk = max(1, -(-n_records // n_parts))  # ceiling division

    chunk_paths = []
    chunk_idx   = 0
    written     = 0
    out_fh      = None

    with open_file(input_fastq) as fh:
        while True:
            name = fh.readline()
            if not name:
                break
            seq  = fh.readline()
            plus = fh.readline()
            qual = fh.readline()

            if written == 0 or written >= records_per_chunk:
                if out_fh:
                    out_fh.close()
                path = os.path.join(tmp_dir, f'chunk_{chunk_idx:04d}.fq.gz')
                chunk_paths.append(path)
                out_fh = gzip.open(path, 'wt')
                chunk_idx += 1
                written = 0

            out_fh.write(name)
            out_fh.write(seq)
            out_fh.write(plus)
            out_fh.write(qual)
            written += 1

    if out_fh:
        out_fh.close()

    return [p for p in chunk_paths if os.path.getsize(p) > 30]


def _process_chunk(records: List[Tuple[str, str, str]]) -> Tuple[List[str], Dict]:
    """Process a list of (read_name, seq, qual) records using per-worker globals."""
    output_lines = []
    stats = {
        'total_reads': 0,
        'reads_with_correctable_bc': 0,
        'reads_with_remapped_bc': 0,
        'mismatch_distribution': defaultdict(int),
    }

    for read_name, seq_line, qual_line in records:
        stats['total_reads'] += 1

        seq_no_sep  = seq_line[1:]  if seq_line.startswith('N')  else seq_line
        qual_no_sep = qual_line[1:] if len(qual_line) == len(seq_line) and seq_line.startswith('N') else qual_line

        right_16_bases = seq_no_sep[-16:]
        right_16_qual  = qual_no_sep[-16:]

        bc_for_matching = right_16_bases[::-1].translate(_RC_TABLE).upper()
        bc_qual = right_16_qual[::-1]  # aligned with bc_for_matching position-wise

        matched_bc, num_mismatches = _match_barcode(bc_for_matching, _wl, _max_mm, bc_qual)

        cr_tag = bc_for_matching
        cy_tag = bc_qual

        if matched_bc is not None:
            stats['reads_with_correctable_bc'] += 1
            stats['mismatch_distribution'][num_mismatches] += 1
            remapped = _remap.get(matched_bc, matched_bc)
            if remapped != matched_bc:
                stats['reads_with_remapped_bc'] += 1
            output_lines.append(f"{read_name} CR:Z:{cr_tag}\tCY:Z:{cy_tag}\tCB:Z:{remapped}-{_suffix}\n")
        else:
            output_lines.append(f"{read_name} CR:Z:{cr_tag}\tCY:Z:{cy_tag}\n")

    return output_lines, stats


def _worker_process_chunk_file(args: Tuple[str, str]) -> Dict:
    """
    Worker entry point for parallel mode.
    Reads chunk_path (gzipped FASTQ), processes every record using the
    per-worker globals loaded by _worker_init, writes text to out_path.
    Returns aggregated stats dict (no record data crosses process boundaries).
    """
    chunk_path, out_path = args
    stats = {
        'total_reads': 0,
        'reads_with_correctable_bc': 0,
        'reads_with_remapped_bc': 0,
        'mismatch_distribution': defaultdict(int),
    }
    with open_file(chunk_path) as fh, open(out_path, 'w') as out_fh:
        for chunk in _read_chunks(fh, 50_000):
            lines, chunk_stats = _process_chunk(chunk)
            out_fh.writelines(lines)
            stats['total_reads']               += chunk_stats['total_reads']
            stats['reads_with_correctable_bc'] += chunk_stats['reads_with_correctable_bc']
            stats['reads_with_remapped_bc']    += chunk_stats['reads_with_remapped_bc']
            for mm, cnt in chunk_stats['mismatch_distribution'].items():
                stats['mismatch_distribution'][mm] += cnt
    return stats


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
        return sequence[::-1].translate(_RC_TABLE).upper()

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
        num_threads: int = 4,
        chunk_size: int = 50_000,
    ) -> None:
        """
        num_threads == 1: single-threaded streaming, zero overhead.
        num_threads  > 1:
          1. Split input FASTQ round-robin into num_threads gzipped chunk
             files in a .bc_tmp/ dir next to the output file.
          2. Run one worker process per chunk — each reads its own file from
             disk and writes its own output file. Zero per-task IPC overhead.
          3. Concatenate output files in order, then delete .bc_tmp/.
        """
        # Set globals for the single-threaded path (workers do it via initializer)
        _worker_init(self.whitelist_barcodes, self.remapping_dict,
                     self.max_mismatches, self.barcode_suffix)

        def _accumulate(s: Dict) -> None:
            self.total_reads               += s['total_reads']
            self.reads_with_correctable_bc += s['reads_with_correctable_bc']
            self.reads_with_remapped_bc    += s['reads_with_remapped_bc']
            for mm, cnt in s['mismatch_distribution'].items():
                self.mismatch_distribution[mm] += cnt

        if num_threads == 1:
            with open_file(input_fastq) as fh:
                for chunk in _read_chunks(fh, chunk_size):
                    lines, stats = _process_chunk(chunk)
                    output_file.writelines(lines)
                    _accumulate(stats)
            return

        # Derive tmp dir from the output file path
        out_name = getattr(output_file, 'name', None)
        if out_name and out_name not in ('<stdout>', '<stderr>'):
            tmp_dir = os.path.join(os.path.dirname(os.path.abspath(out_name)), '.bc_tmp')
        else:
            tmp_dir = os.path.join(os.path.dirname(os.path.abspath(input_fastq)), '.bc_tmp')

        os.makedirs(tmp_dir, exist_ok=True)
        try:
            # Step 1: split
            chunk_paths = _split_fastq(input_fastq, num_threads, tmp_dir)
            out_paths   = [os.path.join(tmp_dir, f'out_{i:04d}.txt') for i in range(len(chunk_paths))]

            # Step 2: process chunks in parallel — workers only pass back small stats dicts
            with Pool(
                processes=len(chunk_paths),
                initializer=_worker_init,
                initargs=(self.whitelist_barcodes, self.remapping_dict,
                          self.max_mismatches, self.barcode_suffix),
            ) as pool:
                all_stats = pool.map(_worker_process_chunk_file, zip(chunk_paths, out_paths))

            # Step 3: concatenate output files in order and accumulate stats
            for out_path, stats in zip(out_paths, all_stats):
                with open(out_path) as fh:
                    shutil.copyfileobj(fh, output_file)
                _accumulate(stats)

        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)

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
                fh.write(f"mismatches_{mismatches}\t{stats['mismatch_distribution'][mismatches]}\n")


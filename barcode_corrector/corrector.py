import gzip
from pathlib import Path
from typing import TextIO, Dict, Tuple, List
from collections import defaultdict
from multiprocessing import Pool
from functools import partial

from .loader import open_file, load_barcode_whitelist, load_barcode_remapping
from .matcher import find_closest_barcode, remap_barcode


def _process_record_batch(batch_records: List[Tuple], 
                          whitelist_barcodes: set, 
                          remapping_dict: dict, 
                          max_mismatches: int,
                          barcode_suffix: str) -> Tuple[List[str], Dict]:
    
    output_lines = []
    stats = {
        'total_reads': 0,
        'reads_with_correctable_bc': 0,
        'reads_with_remapped_bc': 0,
        'mismatch_distribution': defaultdict(int)
    }
    
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    
    for read_name, sequence_line, quality_line in batch_records:
        raw_bc = sequence_line
        raw_bc_qual = quality_line
        stats['total_reads'] += 1
        
        seq_no_sep = raw_bc[1:] if raw_bc.startswith('N') else raw_bc
        qual_no_sep = raw_bc_qual[1:] if len(raw_bc_qual) == len(raw_bc) and raw_bc.startswith('N') else raw_bc_qual
        
        right_16_bases = seq_no_sep[-16:] if len(seq_no_sep) >= 16 else seq_no_sep
        right_16_qual = qual_no_sep[-16:] if len(qual_no_sep) >= 16 else qual_no_sep
        
        bc_for_matching = ''.join(complement.get(base, 'N') for base in reversed(right_16_bases.upper()))
        
        closest_barcode = None
        min_distance = max_mismatches + 1
        best_quality_score = -1
        
        for whitelist_bc in whitelist_barcodes:
            distance = sum(c1 != c2 for c1, c2 in zip(bc_for_matching, whitelist_bc))
            if distance <= max_mismatches:
                if len(right_16_qual) == len(bc_for_matching):
                    quality_score = sum(ord(right_16_qual[i]) for i, (raw_base, ref_base) in enumerate(zip(bc_for_matching, whitelist_bc)) if raw_base != ref_base)
                    if distance < min_distance or (distance == min_distance and quality_score < best_quality_score):
                        closest_barcode = whitelist_bc
                        min_distance = distance
                        best_quality_score = quality_score
                else:
                    if distance < min_distance:
                        closest_barcode = whitelist_bc
                        min_distance = distance
        
        cr_tag = bc_for_matching
        cy_tag = right_16_qual[::-1]
        
        if closest_barcode is not None:
            stats['reads_with_correctable_bc'] += 1
            stats['mismatch_distribution'][min_distance] += 1
            
            remapped_barcode = remapping_dict.get(closest_barcode, closest_barcode)
            if remapped_barcode != closest_barcode:
                stats['reads_with_remapped_bc'] += 1
            
            corrected_bc_with_suffix = f"{remapped_barcode}-{barcode_suffix}"
            output_lines.append(f"{read_name} CR:Z:{cr_tag}\tCY:Z:{cy_tag}\tCB:Z:{corrected_bc_with_suffix}\n")
        else:
            output_lines.append(f"{read_name} CR:Z:{cr_tag}\tCY:Z:{cy_tag}\n")
    
    return output_lines, stats


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
        num_threads: int = 4
    ) -> None:
        
        fastq_records = []
        with open_file(input_fastq) as fh:
            while True:
                read_name_line = fh.readline().strip()
                if not read_name_line:
                    break
                
                sequence_line = fh.readline().strip()
                plus_line = fh.readline().strip()
                quality_line = fh.readline().strip()
                
                read_name = read_name_line[1:] if read_name_line.startswith('@') else read_name_line
                fastq_records.append((read_name, sequence_line, quality_line))
        
        if not fastq_records:
            return
        
        chunk_size = max(1, len(fastq_records) // num_threads)
        chunks = [fastq_records[i:i + chunk_size] for i in range(0, len(fastq_records), chunk_size)]
        
        process_func = partial(
            _process_record_batch,
            whitelist_barcodes=self.whitelist_barcodes,
            remapping_dict=self.remapping_dict,
            max_mismatches=self.max_mismatches,
            barcode_suffix=self.barcode_suffix
        )
        
        with Pool(processes=num_threads) as pool:
            results = pool.map(process_func, chunks)
        
        for output_lines, chunk_stats in results:
            for line in output_lines:
                output_file.write(line)
            
            self.total_reads += chunk_stats['total_reads']
            self.reads_with_correctable_bc += chunk_stats['reads_with_correctable_bc']
            self.reads_with_remapped_bc += chunk_stats['reads_with_remapped_bc']
            for mismatches, count in chunk_stats['mismatch_distribution'].items():
                self.mismatch_distribution[mismatches] += count

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

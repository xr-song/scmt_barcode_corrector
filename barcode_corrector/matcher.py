from typing import Optional, Set


def hamming_distance(seq1: str, seq2: str) -> int:
    if len(seq1) != len(seq2):
        return max(len(seq1), len(seq2))
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def find_closest_barcode(
    raw_barcode: str,
    whitelist_barcodes: Set[str],
    max_mismatches: int,
    quality: Optional[str] = None
) -> Optional[str]:
    closest_barcode = None
    min_distance = max_mismatches + 1
    best_quality_score = -1

    for whitelist_bc in whitelist_barcodes:
        distance = hamming_distance(raw_barcode, whitelist_bc)
        if distance <= max_mismatches:
            if quality and len(quality) == len(raw_barcode):
                quality_score = 0
                for i, (raw_base, ref_base) in enumerate(zip(raw_barcode, whitelist_bc)):
                    if raw_base != ref_base:
                        quality_score += ord(quality[i])
                
                if distance < min_distance or (distance == min_distance and quality_score < best_quality_score):
                    closest_barcode = whitelist_bc
                    min_distance = distance
                    best_quality_score = quality_score
            else:
                if distance < min_distance:
                    closest_barcode = whitelist_bc
                    min_distance = distance

    return closest_barcode


def remap_barcode(
    corrected_barcode: str,
    remapping_dict: dict
) -> str:
    return remapping_dict.get(corrected_barcode, corrected_barcode)

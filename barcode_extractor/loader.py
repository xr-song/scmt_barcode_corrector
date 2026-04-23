import gzip
from pathlib import Path
from typing import Set, Dict, TextIO


def open_file(filename: str) -> TextIO:
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r')


def load_barcode_whitelist(whitelist_file: str) -> Set[str]:
    barcodes = set()
    with open_file(whitelist_file) as fh:
        for line in fh:
            line = line.strip()
            if line:
                barcode = line.split('\t')[0]
                barcodes.add(barcode)
    return barcodes


def load_barcode_remapping(whitelist_file: str, remapping_file: str) -> Dict[str, str]:
    if not remapping_file or remapping_file.lower() in ('false', 'none'):
        return {}

    remapping = {}
    with open_file(whitelist_file) as wh, open_file(remapping_file) as rh:
        for whitelist_line, remapping_line in zip(wh, rh):
            whitelist_bc = whitelist_line.strip().split('\t')[0]
            remapping_bc = remapping_line.strip().split('\t')[0]
            if whitelist_bc and remapping_bc:
                remapping[whitelist_bc] = remapping_bc
    return remapping


def validate_barcode_whitelist(whitelist_file: str) -> int:
    with open_file(whitelist_file) as fh:
        first_line = fh.readline().strip()
        if not first_line:
            raise ValueError(f"Whitelist file {whitelist_file} is empty")

        first_barcode = first_line.split('\t')[0]

        if not first_barcode:
            raise ValueError(f"No valid barcode in first line of {whitelist_file}")

        if not all(c in 'ACGT' for c in first_barcode):
            raise ValueError(
                f"First barcode '{first_barcode}' contains non-ACGT characters"
            )

        return len(first_barcode)

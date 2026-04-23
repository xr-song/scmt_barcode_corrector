# barcode-extractor

Efficient Python implementation of barcode correction for in-house single-cell genomics FASTQ files.

## Barcode correction logic

The barcode corrector processes FASTQ files containing raw barcode sequences from single-cell genomics experiments:

```mermaid
flowchart TD
    A["read fastq record"] --> B["extract sequence\n& quality"]
    B --> D["extract last 16 bases\n& reverse complement"]
    D --> F["find closest match in DNA whitelist"]
    F --> G{"match found?"}
    G -->|yes| H["map to RNA whitelist"]
    G -->|no| I["keep raw barcode"]
    H --> J["output: CR:Z, CY:Z, CB:Z tags"]
    I --> K["output: CR:Z, CY:Z tags"]
    J --> L["write barcode detection stat"]
    K --> L
```

## Installation

```bash
git clone <this repo>
pip install -e .
```

## Usage

```bash
correct-barcodes <dna_whitelist_file> <rna_whitelist_file> <input_fastq> <output_file> [options]
```

### Example

```bash
correct-barcodes \
    737K-arc-v1.txt.gz \
    737K-rna-v1.txt.gz \
    raw_barcodes.fastq.gz \
    corrected_barcodes.tsv \
    --barcode-suffix 1 \
    --max-mismatches 1 \
    --stats-file correction_stats.tsv \
    --threads 8
```

## Parameters

- `dna_whitelist_file`: DNA barcode whitelist (one per line)
- `rna_whitelist_file`: RNA barcode whitelist for remapping. Use "false" or "none" to disable.
- `input_fastq`: Input FASTQ with raw barcodes (gzip or uncompressed)
- `output_file`: Output TSV with corrected barcodes
- `--barcode-suffix`: Suffix to append (default: 1)
- `--max-mismatches`: Max mismatches (default: 1)
- `--min-fraction`: Min correctable fraction (default: 0.5)
- `--stats-file`: Statistics output file
- `--threads`: CPUs to use (default: 1)

## Output format

Tab-separated columns:
1. read_name
2. CR:Z - raw barcode sequence (reverse complemented)
3. CY:Z - raw barcode quality (reversed to match rc)
4. CB:Z - corrected barcode with suffix (if correctable)

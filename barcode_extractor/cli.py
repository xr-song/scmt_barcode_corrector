import argparse
import gzip
import sys
from pathlib import Path

from .loader import validate_barcode_whitelist
from .corrector import BarcodeCorrector


def main():
    parser = argparse.ArgumentParser(
        description="Correct raw barcodes from FASTQ index files against a whitelist."
    )

    parser.add_argument(
        'dna_whitelist_file',
        help='DNA barcode whitelist file'
    )

    parser.add_argument(
        'rna_whitelist_file',
        help='RNA barcode whitelist file'
    )

    parser.add_argument(
        'input_fastq',
        help='Input FASTQ index file with raw barcodes (gzipped or uncompressed)'
    )

    parser.add_argument(
        'output_file',
        help='Output file with corrected barcodes (supports .gz and .zst compression)'
    )

    parser.add_argument(
        '--barcode-suffix',
        default='1',
        help='Suffix to append after corrected barcode (default: 1)'
    )

    parser.add_argument(
        '--max-mismatches',
        type=int,
        default=1,
        help='Maximum mismatches allowed for barcode correction (default: 1)'
    )

    parser.add_argument(
        '--min-fraction',
        type=float,
        default=0.5,
        help='Minimum fraction of reads with correctable barcodes (default: 0.5)'
    )

    parser.add_argument(
        '--stats-file',
        help='Output file for correction statistics (optional)'
    )

    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of threads to use for compression (default: 1)'
    )

    args = parser.parse_args()

    try:
        if not Path(args.dna_whitelist_file).exists():
            print(
                f"Error: DNA whitelist file '{args.dna_whitelist_file}' not found",
                file=sys.stderr
            )
            return 1

        if not Path(args.input_fastq).exists():
            print(
                f"Error: Input FASTQ file '{args.input_fastq}' not found",
                file=sys.stderr
            )
            return 1

        barcode_length = validate_barcode_whitelist(args.dna_whitelist_file)
        print(f"Loaded barcode whitelist (barcode length: {barcode_length} bp)")

        corrector = BarcodeCorrector(
            dna_whitelist_file=args.dna_whitelist_file,
            rna_whitelist_file=args.rna_whitelist_file,
            max_mismatches=args.max_mismatches,
            min_frac_bcs_to_find=args.min_fraction,
            barcode_suffix=args.barcode_suffix
        )

        print("Starting barcode extraction...")

        output_path = args.output_file
        use_gzip = output_path.endswith('.gz')
        use_zst = output_path.endswith('.zst')
        
        if use_zst:
            output_path_tmp = output_path[:-4]
        elif use_gzip:
            output_path_tmp = output_path[:-3]
        else:
            output_path_tmp = output_path

        with open(output_path_tmp, 'w') as out_fh:
            corrector.process_fastq(args.input_fastq, out_fh, num_threads=args.threads)

        if use_gzip:
            with open(output_path_tmp, 'rb') as f_in:
                with gzip.open(args.output_file, 'wb') as f_out:
                    f_out.write(f_in.read())
            Path(output_path_tmp).unlink()
        
        if use_zst:
            try:
                import zstandard as zstd
                with open(output_path_tmp, 'rb') as f_in:
                    with open(args.output_file, 'wb') as f_out:
                        cctx = zstd.ZstdCompressor(level=6, threads=args.threads)
                        f_out.write(cctx.compress(f_in.read()))
                Path(output_path_tmp).unlink()
            except ImportError:
                print(
                    "Warning: zstandard not installed. Keeping uncompressed file.",
                    file=sys.stderr
                )

        stats = corrector.get_statistics()
        print(f"\nCorrection Statistics:")
        print(f"  Total reads: {stats['total_reads']}")
        print(f"  Reads with correctable barcode: {stats['reads_with_correctable_bc']}")
        print(f"  Fraction correctable: {stats['fraction_correctable']:.2%}")
        print(f"  Reads with remapped barcode: {stats['reads_with_remapped_bc']}")

        if args.stats_file:
            corrector.write_statistics(args.stats_file)
            print(f"  Statistics written to: {args.stats_file}")

        print("\nBarcode correction complete!")
        return 0

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

import argparse
from main_dna_rna_filter import FastqFilter


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter FASTQ records based on GC content, length, and quality."
        )
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to input FASTQ file"
        )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to output FASTQ file"
        )
    parser.add_argument(
        "--min_gc",
        type=float,
        default=0.0,
        help="Minimum GC content (default: 0)"
                        )
    parser.add_argument(
        "--max_gc",
        type=float,
        default=100.0,
        help="Maximum GC content (default: 100)"
        )
    parser.add_argument(
        "--min_len",
        type=int,
        default=0,
        help="Minimum read length (default: 0)"
        )
    parser.add_argument(
        "--max_len",
        type=int,
        default=2**32,
        help="Maximum read length (default: 2^32)"
        )
    parser.add_argument(
        "--min_quality",
        type=float,
        default=0.0,
        help="Minimum average quality (default: 0)"
        )
    return parser.parse_args()


def main():
    args = parse_args()

    filter_tool = FastqFilter(
        gc_bounds=(args.min_gc, args.max_gc),
        length_bounds=(args.min_len, args.max_len),
        quality_threshold=args.min_quality
    )

    filter_tool.filter_fastq(args.input, args.output)


if __name__ == "__main__":
    main()

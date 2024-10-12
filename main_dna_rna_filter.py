from utils.dna_rna_helper import (
    transcribe, reverse, complement, reverse_complement, gc_content,
    is_palindrome, melting_temperature, change_sequence
)
from utils.filter_helper import avg_quality, read_fastq as fastq_to_dict, \
    write_fastq, gc_content as filter_gc_content
from typing import Tuple, Any, List, Union


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Tuple[float, float] = (0, 100),
    length_bounds: Tuple[int, int] = (0, 2**32),
    quality_threshold: float = 0
) -> None:
    """
    Filters reads based on GC content, sequence length, and average quality,
    and writes filtered reads to a file.

    Arguments:
    - input_fastq: Path to the input FASTQ file.
    - output_fastq: Path to the output FASTQ file.
    - gc_bounds: Range for the percentage of GC content (default is (0, 100)).
    - length_bounds: Range for the sequence length (default is (0, 2**32)).
    - quality_threshold: Threshold for the average quality of the read
      (default is 0).

    Returns:
    - None. Writes the filtered sequences to the output FASTQ file.
    """

    seqs = fastq_to_dict(input_fastq)

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    filtered_seqs = {}

    for name, (seq, quality) in seqs.items():
        gc = filter_gc_content(seq)
        length = len(seq)
        avg_qual = avg_quality(quality)

        if (gc_bounds[0] <= gc <= gc_bounds[1] and
                length_bounds[0] <= length <= length_bounds[1] and
                avg_qual >= quality_threshold):
            filtered_seqs[name] = (seq, quality)

    write_fastq(output_fastq, filtered_seqs)


def run_dna_rna_tools(*args: str, **kwargs: Any) -> Union[str, List[str]]:
    """
    Performs various DNA/RNA sequence operations such as transcription,
    complement, reverse complement, and more.

    Arguments:
    - *args: The DNA/RNA sequences and the operation to perform.
      The last argument must be the operation name (str).
        - Example operations: "transcribe", "complement", "reverse_complement",
          "gc_content", "is_palindrome", "melting_temperature",
          "change_sequence", "reverse".
    - **kwargs: Additional keyword arguments to customize operations
      (e.g., "case" for case transformation in `change_sequence`).

    Returns:
    - A string result of the operation for a single sequence or a list of
      results for multiple sequences. If the sequences contain invalid
      characters, an error message will be returned.
    """

    def validate_seq(seq: str) -> str:
        if ("T" in seq or "t" in seq) and ("U" in seq or "u" in seq):
            return "Error: Sequence contains both T and U."

        if "T" in seq or "t" in seq:
            valid_set = set("ATGCatgc")
        else:
            valid_set = set("AUGCaugc")

        if set(seq) - valid_set:
            return "Error: Sequence contains invalid characters."

        return None

    if len(args) < 2:
        return "Invalid number of arguments."

    *sequences, operation = args

    operations = {
        "transcribe": transcribe,
        "complement": complement,
        "reverse_complement": reverse_complement,
        "gc_content": gc_content,
        "is_palindrome": is_palindrome,
        "melting_temperature": melting_temperature,
        "change_sequence": change_sequence,
        "reverse": reverse,
    }

    if operation not in operations:
        return "Operation not supported."

    results = []
    for seq in sequences:
        validation_error = validate_seq(seq)
        if validation_error:
            return validation_error
        result = operations[operation](seq, **kwargs)
        results.append(result)

    if len(results) == 1:
        return results[0]
    return results

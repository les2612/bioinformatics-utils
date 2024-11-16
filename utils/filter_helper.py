import os
from typing import Dict, Tuple


def gc_content(seq: str) -> float:
    """
    Calculates the percentage of GC content in a sequence as a float.
    Returns 0.0 if the sequence is empty.
    """
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq) * 100) if seq else 0.0


def avg_quality(quality_str: str) -> float:
    """
    Calculates the average quality score of a sequence based on ASCII codes
    (Phred33 scale) as a float. Returns 0.0 if the quality string is empty.
    """
    if not quality_str:
        return 0.0
    return sum(ord(char) - 33 for char in quality_str) / len(quality_str)


def read_fastq(input_fastq: str) -> Dict[str, Tuple[str, str]]:
    """
    Reads a FASTQ file and converts it into a dictionary.

    Arguments:
    - input_fastq: Path to the input FASTQ file.

    Returns:
    - A dictionary where the key is the sequence header, and the value is a
      tuple (sequence, quality).
    """
    fastq_dict = {}
    with open(input_fastq, 'r') as file:
        line_count = 0
        current_entry = []

        for line in file:
            line = line.strip()
            current_entry.append(line)
            line_count += 1

            if line_count == 4:
                name = current_entry[0]
                sequence = current_entry[1]
                quality = current_entry[3]
                fastq_dict[name] = (sequence, quality)
                current_entry = []
                line_count = 0

    return fastq_dict


def write_fastq(output_fastq: str, seqs: Dict[str, Tuple[str, str]]):
    """
    Writes filtered sequences to the output FASTQ file in the 'filtered'
    folder.

    Arguments:
    - output_fastq: Name of the output FASTQ file.
    - seqs: Dictionary where the key is the sequence header, and the value is
      a tuple (sequence, quality).

    The function creates the 'filtered' folder if it does not exist.
    """
    output_dir = 'filtered'
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, output_fastq)

    with open(output_path, 'w') as file:
        for name, (sequence, quality) in seqs.items():
            file.write(name + "\n")
            file.write(sequence + "\n")
            file.write("+\n")
            file.write(quality + "\n")

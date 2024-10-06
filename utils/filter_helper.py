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

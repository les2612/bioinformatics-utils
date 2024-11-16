from typing import Any


def transcribe(seq: str) -> str:
    """
    Converts DNA sequence to RNA by replacing T with U.
    """
    return seq.replace('T', 'U').replace('t', 'u')


def reverse(seq: str) -> str:
    """
    Reverses the given sequence.
    """
    return seq[::-1]


def complement(seq: str) -> str:
    """
    Returns the complement of the DNA or RNA sequence.
    """
    dna_rule = {
        "A": "T", "T": "A", "G": "C", "C": "G",
        "a": "t", "t": "a", "g": "c", "c": "g",
    }

    rna_rule = {
        "A": "U", "U": "A", "G": "C", "C": "G",
        "a": "u", "u": "a", "g": "c", "c": "g",
    }

    rule = rna_rule if "U" in seq or "u" in seq else dna_rule

    result = [rule.get(chair, chair) for chair in seq]
    return "".join(result)


def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complement of a DNA or RNA sequence.
    """
    return reverse(complement(seq))


def gc_content(seq: str) -> float:
    """
    Calculates the GC content (percentage) of a DNA or RNA sequence.
    """
    gc_count = 0
    for chair in seq:
        if chair.upper() == "G" or chair.upper() == "C":
            gc_count += 1
    return (gc_count / len(seq)) * 100 if seq else 0.0


def is_palindrome(seq: str) -> bool:
    """
    Checks if the sequence is a palindrome.
    """
    return seq == reverse(seq)


def melting_temperature(seq: str) -> int:
    """
    Calculates the melting temperature of a DNA sequence
    using the basic formula.
    """
    a_count = seq.count("A") + seq.count("a")
    t_count = seq.count("T") + seq.count("t")
    g_count = seq.count("G") + seq.count("g")
    c_count = seq.count("C") + seq.count("c")

    return 4 * (g_count + c_count) + 2 * (a_count + t_count)


def change_sequence(seq: str, **kwargs: Any) -> str:
    """
    Changes the case of the sequence based on kwargs (uppercase or lowcase).
    """
    if "case" in kwargs and kwargs["case"] == "uppercase":
        return seq.upper()
    if "case" in kwargs and kwargs["case"] == "lowcase":
        return seq.lower()
    return seq

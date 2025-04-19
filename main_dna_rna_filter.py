from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord
from typing import Tuple, List
import logging
import os


os.makedirs("logs", exist_ok=True)
logging.basicConfig(
    filename="logs/filter.log",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)


class BiologicalSequence(ABC):
    """Abstract base class for all biological sequences."""

    def __init__(self, sequence: str):
        """
        Initialize a biological sequence.

        :param sequence: Nucleotide or amino acid sequence string.
        """
        self.sequence = sequence

    @abstractmethod
    def __len__(self) -> int:
        """
        Return the length of the sequence.
        """
        pass

    @abstractmethod
    def __getitem__(self, index):
        """
        Return the nucleotide or amino acid at the specified index.
        """
        pass

    @abstractmethod
    def __str__(self) -> str:
        """
        Return a string representation of the sequence.
        """
        pass

    @abstractmethod
    def check_alphabet(self) -> bool:
        """
        Check whether the sequence contains only valid characters.
        """
        pass


class MolecularSequence(BiologicalSequence):
    """Intermediate class for molecular sequences (DNA, RNA, proteins)."""
    valid_characters = set()

    def check_alphabet(self) -> bool:
        """
        Check if all characters in the sequence belong to the valid
        character set.
        """
        return set(self.sequence.upper()).issubset(self.valid_characters)


class NucleicAcidSequence(MolecularSequence):
    """Represents a nucleic acid sequence (DNA or RNA)."""
    complement_map = {}

    def __init__(self, sequence: str):
        """
    Initialize a nucleic acid sequence.

    This constructor prevents direct instantiation of the abstract base class
    NucleicAcidSequence. It should only be called from a subclass
    (e.g. DNASequence or RNASequence).

    :param sequence: The nucleotide sequence string.
    :raises NotImplementedError: If attempting to instantiate
    NucleicAcidSequence directly.
    """
        if self.__class__ == NucleicAcidSequence:
            raise NotImplementedError(
                "NucleicAcidSequence is an abstract class and cannot be "
                "instantiated directly."
            )
        super().__init__(sequence)

    def complement(self) -> "NucleicAcidSequence":
        """
        Return the complement of the sequence.

        :return: Complemented sequence as the same subclass.
        """
        return self.__class__(
            "".join(self.complement_map.get(base, base)
                    for base in self.sequence)
            )

    def reverse(self) -> "NucleicAcidSequence":
        """
        Return the reversed sequence.
        """
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        """
        Return the reverse complement of the sequence.
        """
        return self.complement().reverse()

    def gc_content(self) -> float:
        """
        Calculate the GC content (%) of the sequence.

        :return: GC content as a float percentage.
        """
        if not self.sequence:
            raise ValueError(
                "GC content cannot be calculated for an empty sequence."
                )
        return gc_fraction(self.sequence) * 100

    def melting_temperature(self) -> int:
        """
        Estimate the melting temperature of the sequence using the
        Wallace rule.

        :return: Melting temperature in Celsius.
        """
        if not self.sequence:
            raise ValueError(
                "Melting temperature cannot be calculated "
                "for an empty sequence."
                )
        return (
            4 * (self.sequence.count("G") + self.sequence.count("C")) +
            2 * (self.sequence.count("A") + self.sequence.count("T"))
        )


class DNASequence(NucleicAcidSequence):
    """Represents a DNA sequence."""
    complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    valid_characters = {"A", "T", "G", "C"}

    def transcribe(self) -> "RNASequence":
        """
        Transcribe the DNA sequence into RNA.

        :return: Corresponding RNASequence with 'T' replaced by 'U'.
        """
        return RNASequence(self.sequence.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    """Represents an RNA sequence."""
    complement_map = {"A": "U", "U": "A", "G": "C", "C": "G"}
    valid_characters = {"A", "U", "G", "C"}


class AminoAcidSequence(MolecularSequence):
    """Represents an amino acid sequence (protein)."""
    valid_characters = set("ACDEFGHIKLMNPQRSTVWY")

    def molecular_weight(self) -> float:
        """
        Calculate the approximate molecular weight of the protein sequence.

        :return: Molecular weight in Daltons (Da).
        """
        amino_acid_weights = {
            "A": 89.1, "C": 121.2, "D": 133.1, "E": 147.1,
            "F": 165.2, "G": 75.1, "H": 155.2, "I": 131.2,
            "K": 146.2, "L": 131.2, "M": 149.2, "N": 132.1,
            "P": 115.1, "Q": 146.2, "R": 174.2, "S": 105.1,
            "T": 119.1, "V": 117.1, "W": 204.2, "Y": 181.2
        }
        return sum(
            amino_acid_weights.get(aa, 0) for aa in self.sequence.upper()
            )


class FastqFilter:
    """Class for filtering FASTQ files based on GC content, length,
    and quality."""

    def __init__(
        self,
        gc_bounds: Tuple[float, float] = (0, 100),
        length_bounds: Tuple[int, int] = (0, 2**32),
        quality_threshold: float = 0,
    ):
        """
        Initializes the filter with specified thresholds.

        :param gc_bounds: GC content range (percentage).
        :param length_bounds: Sequence length range.
        :param quality_threshold: Minimum average quality score.
        """
        if gc_bounds[0] > gc_bounds[1]:
            raise ValueError(
                f"Invalid GC bounds: {gc_bounds}."
                f"Lower bound cannot be greater than upper bound."
                )
        self.gc_bounds = gc_bounds
        self.length_bounds = length_bounds
        self.quality_threshold = quality_threshold

    @staticmethod
    def avg_quality(quality_scores: List[int]) -> float:
        """
        Computes the average quality score of a sequence.

        :param quality_scores: List of Phred quality scores.
        :return: Average quality score.
        """
        if not quality_scores:
            raise ValueError(
                "Cannot compute average quality: "
                "no quality scores provided.")
        return sum(quality_scores) / len(quality_scores)

    def filter_records(self, records: List[SeqRecord]) -> List[SeqRecord]:
        """
        Filters records based on GC content, length, and quality.

        :param records: List of SeqRecord objects.
        :return: List of filtered SeqRecord objects.
        """
        filtered_records = []
        for record in records:
            try:
                if not record.seq or len(record.seq) == 0:
                    raise ValueError(f"Empty sequence found in "
                                     f"record: {record.id}")
                gc_content = gc_fraction(record.seq) * 100
                if not (
                    self.gc_bounds[0] <= gc_content <= self.gc_bounds[1]
                ):
                    raise ValueError(
                        f"GC content {gc_content:.2f}% out of bounds "
                        f"for record: {record.id}"
                    )

                if not (
                    self.length_bounds[0] <= len(record.seq)
                    <= self.length_bounds[1]
                ):
                    raise ValueError(
                        f"Sequence length {len(record.seq)} out of "
                        f"bounds for record: {record.id}"
                    )

                if "phred_quality" not in record.letter_annotations:
                    raise ValueError(f"Missing Phred quality scores "
                                     f"in record: {record.id}")

                if not record.letter_annotations["phred_quality"]:
                    raise ValueError(f"Empty Phred quality scores "
                                     f"in record: {record.id}")

                avg_quality = self.avg_quality(
                    record.letter_annotations["phred_quality"]
                    )
                if avg_quality < self.quality_threshold:
                    raise ValueError(
                        f"Average quality {avg_quality:.2f} below threshold "
                        f"for record: {record.id}"
                    )

                filtered_records.append(record)

            except ValueError as e:
                error_message = f"Record {record.id}: {e}"
                print(f"Warning: {error_message}")
                logging.error(error_message)

        return filtered_records

    def filter_fastq(self, input_fastq: str, output_fastq: str):
        """
        Loads, filters, and saves a FASTQ file.

        :param input_fastq: Path to the input FASTQ file.
        :param output_fastq: Path to save the filtered FASTQ file.
        """
        records = list(SeqIO.parse(input_fastq, "fastq"))
        if not records:
            raise ValueError(f"No records found in FASTQ file: {input_fastq}")
        filtered_records = self.filter_records(records)
        self.save_filtered(output_fastq, filtered_records)

    @staticmethod
    def save_filtered(output_fastq: str, filtered_records: List[SeqRecord]):
        """
        Saves filtered sequences to a FASTQ file.

        :param output_fastq: Path to the output file.
        :param filtered_records: List of filtered sequences.
        """
        SeqIO.write(filtered_records, output_fastq, "fastq")
        message = (
            f"Filtered {len(filtered_records)} "
            f"reads saved to {output_fastq}."
            )
        print(message)
        logging.info(message)

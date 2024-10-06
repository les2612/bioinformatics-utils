# bioinformatics-utils

**bioinformatics-utils** is a comprehensive toolkit for performing basic operations and analyses on DNA and RNA sequences. It provides functionality for filtering FASTQ sequences based on GC content, sequence length, and quality score, as well as utilities for sequence transformations like transcription, reverse complement, palindrome checking, and more.

## Authors:
- **Software Development**: [Your Name], [Your Institution]
- **Testing**: [Collaborator Name], [Institution]
- **Supervisor**: [Supervisor Name], [Institution]

---


## Content
- [Installation](#installation)
- [Running Instructions](#running-instructions)
- [Examples](#examples)
- [FAQ](#faq)
- [Citation](#citation)
- [Contact](#contact)

---

## Installation

To install and use **bioinformatics-utils**, follow the steps below:

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/bioinformatics-utils.git
   cd bioinformatics-utils
   ```

2. Make sure Python 3.8+ is installed, along with pip:
    ```bash
    pip install -r requirements.txt
    ```

3. Now, you can import and use bioinformatics-utils in your Python scripts or Jupyter notebooks.

Alternatively, you can use this package as a standalone script with basic command-line options.

## Running Instructions
**bioinformatics-utils** offers two main functionalities:

1. **FASTQ Filtering:**
The FASTQ Filtering functionality allows you to filter reads based on the following criteria:

*GC Content:* The percentage of guanine and cytosine (GC) in the sequence.
*Sequence Length:* The length of the nucleotide sequence.
*Quality Score:* The average Phred33 quality score of the sequence.

**How to Run:**
You can run the filtering function filter_fastq() by providing it with a dictionary of sequences, and optionally, you can specify filtering parameters such as gc_bounds, length_bounds, and quality_threshold

**Available Parameters:**
gc_bounds (tuple): The minimum and maximum percentage of GC content to filter by (e.g., (40, 60)).
length_bounds (tuple): The minimum and maximum length of sequences to include (e.g., (50, 100)).
quality_threshold (float): The minimum average quality score for a sequence to be included (default is 0).

2. **DNA/RNA Sequence Manipulation:**

The DNA/RNA Sequence Manipulation functionality provides several operations to transform and analyze DNA or RNA sequences:

*Transcription:* Converts a DNA sequence to an RNA sequence (replacing T with U).
*Reverse Complement:* Computes the reverse complement of a sequence.
*Palindrome Detection:* Checks if a sequence is a palindrome.
*Melting Temperature:* Calculates the melting temperature of a DNA sequence based on the nucleotide composition.

**How to Run:**
You can run the run_dna_rna_tools() function with various operations by providing the sequences and specifying the operation as the last argument.

**Available Operations:**
transcribe: Converts DNA to RNA.
complement: Returns the complement of a DNA or RNA sequence.
reverse_complement: Returns the reverse complement of a DNA or RNA sequence.
gc_content: Calculates the percentage of GC content in the sequence.
is_palindrome: Checks if the sequence is a palindrome.
melting_temperature: Calculates the melting temperature of a DNA sequence.
change_sequence: Changes the case of the sequence (e.g., uppercase, lowcase).
reverse: Reverses the sequence.

## Examples

1. FASTQ Filtering:

```python
from bioinformatics_utils import filter_fastq

sequences = {
    "@SEQ_ID1": ("ATGCGTACGTAGCT", "IIIIIIIIIIIIII"),
    "@SEQ_ID2": ("GGCTTACGATCGA", "HHHHHHHHHHHHHH"),
}


filtered = filter_fastq(sequences, gc_bounds=(40, 60), quality_threshold=30)
print(filtered)
```

2. DNA/RNA Sequence Manipulation:

```python
from bioinformatics_utils import filter_fastq

sequences = {
    "@SEQ_ID1": ("ATGCGTACGTAGCT", "IIIIIIIIIIIIII"),
    "@SEQ_ID2": ("GGCTTACGATCGA", "HHHHHHHHHHHHHH"),
}

filtered = filter_fastq(sequences, gc_bounds=(40, 60), quality_threshold=30)
print(filtered)
```


## FAQ
Q: Can bioinformatics-utils handle both uppercase and lowercase sequences?
A: Yes, bioinformatics-utils works seamlessly with both uppercase and lowercase nucleotide sequences.

Q: How do I filter sequences by GC content and quality score in one step?
A: You can pass both the gc_bounds and quality_threshold parameters to filter_fastq() to apply both filters at once.

Q: What types of sequences are supported?
A: bioinformatics-utils supports both DNA and RNA sequences for various operations like transcription and complement generation.

## Citation
If you use bioinformatics-utils in your research, please cite this repository and the following publication:

Lebedeva ES, (2024). bioinformatics-utils: A comprehensive toolkit for DNA and RNA sequence analysis.

## Contact

For any issues or questions, please report directly to the GitHub issue tracker, or contact the author via email at Ekaterinalebedeva2612@gmail.com.


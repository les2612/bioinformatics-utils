# bioinformatics-utils

**bioinformatics-utils** bioinformatics-utils is a comprehensive toolkit for performing basic operations and analyses on DNA and RNA sequences. It provides functionality for filtering FASTQ sequences based on GC content, sequence length, and quality score, as well as utilities for sequence transformations like transcription, reverse complement, palindrome checking, and more. It also includes utility functions for processing FASTA and BLAST output files.

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

3. **Convert FASTA Multiline to Single Line:**
The function convert_multiline_fasta_to_oneline() reads a multi-line FASTA file where sequences may be broken into multiple lines and converts it to a single-line FASTA format, where each sequence appears on a single line.

**How to Run:**
You can provide an input FASTA file to convert and optionally provide an output file for saving the result. If the output file is not provided, the result will be printed to the terminal.

4. **Parse BLAST Output:**
The function parse_blast_output extracts the top protein name from each query in BLAST output files.

## Examples

1. FASTQ Filtering:

```python
from bioinformatics_utils import filter_fastq

input_fastq = "input_file.fastq"
output_fastq = "filtered_output.fastq"
filter_fastq(input_fastq, output_fastq, gc_bounds=(40, 60), length_bounds=(50, 200), quality_threshold=30)
```

2. DNA/RNA Sequence Manipulation:

```python
from bioinformatics_utils import run_dna_rna_tools

result = run_dna_rna_tools("ATGCGTACG", "transcribe")
print(result)  

result = run_dna_rna_tools("ATGCGTACG", "reverse_complement")
print(result)  

result = run_dna_rna_tools("ATCGAT", "is_palindrome")
print(result)  
```

3. Convert FASTA Multiline to Single Line:

The bio_files_processor.py is part of the bioinformatics-utils package and provides two key functions: convert_multiline_fasta_to_oneline() 

```python
from bioinformatics_utils.bio_files_processor import convert_multiline_fasta_to_oneline

input_fasta = "example_multiline_fasta.fasta"
output_fasta = "output.fasta"

convert_multiline_fasta_to_oneline(input_fasta, output_fasta)
```

4. Parse BLAST Output:

```python
from bioinformatics_utils.bio_files_processor import parse_blast_output

input_file = "example_blast_results.txt"
output_file = "protein_list.txt"

parse_blast_output(input_file, output_file)
```
5. Running through the terminal
The script bio_files_processor.py contains two functions: convert_multiline_fasta_to_oneline and parse_blast_output, which can be selected and run through the terminal

```bash
python bio_files_processor.py
```
Follow the prompt to select the function you want to run:

- Enter 1 to run convert_multiline_fasta_to_oneline.
- Enter 2 to run parse_blast_output.


## FAQ
Q: Can bioinformatics-utils handle both uppercase and lowercase sequences?
A: Yes, bioinformatics-utils works seamlessly with both uppercase and lowercase nucleotide sequences.

Q: How do I filter sequences by GC content and quality score in one step?
A: You can pass both the gc_bounds and quality_threshold parameters to filter_fastq() to apply both filters at once.

Q: What types of sequences are supported?
A: bioinformatics-utils supports both DNA and RNA sequences for various operations like transcription and complement generation.

Q: Can I convert multi-line FASTA files to a single-line format?
A: Yes, by using the convert_multiline_fasta_to_oneline() function, you can convert multi-line FASTA files so that each sequence appears on a single line. This function can be run through the terminal or used programmatically within your scripts.

Q: How do I extract the top protein names from BLAST output files?
A: You can use the parse_blast_output() function to extract the top protein names for each query from a BLAST output file. The function supports exporting the results to a file or printing them to the terminal.

Q: Is it possible to run bioinformatics-utils functions directly from the terminal?
A: Yes! You can run filter_fastq, convert_multiline_fasta_to_oneline, and parse_blast_output directly from the terminal. We've provided example scripts in the documentation that show how to set up these functions to run from the command line.


## Citation
If you use bioinformatics-utils in your research, please cite this repository and the following publication:

Lebedeva ES, (2024). bioinformatics-utils: A comprehensive toolkit for DNA and RNA sequence analysis.

## Contact

For any issues or questions, please report directly to the GitHub issue tracker, or contact the author via email at Ekaterinalebedeva2612@gmail.com.

# bioinformatics-utils

**bioinformatics-utils** is a toolkit providing essential utilities for processing biological sequence data. This repository includes two main scripts:
1. **bio_files_processor.py**: Functions for processing FASTA and BLAST output files.
2. **main_dna_rna_filter.py**: Utilities for sequence manipulation and FASTQ filtering based on GC content, sequence length, and quality score.


## Authors:
- **Software Development**: Ekaterina Lebedeva, Bioinformatics Institute
- **Testing**: Ekaterina Lebedeva, Bioinformatics Institute
- **Supervisor**: Nikita Vaulin, Bioinformatics Institute

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
Functions:

- convert_multiline_fasta_to_oneline: Converts multi-line FASTA files into single-line format.
- parse_blast_output: Extracts the top protein name from each query in BLAST output files.

Running from Terminal:
```bash
python bio_files_processor.py
```

Follow the prompt to select:

- Enter 1 for FASTA conversion.
- Enter 2 for BLAST output parsing.

2. . Sequence Manipulation & FASTQ Filtering (main_dna_rna_filter.py)

Features:

- Transcription: DNA → RNA conversion.
- Reverse Complement: Computes the reverse complement of a sequence.
- GC Content Calculation: Computes the GC content of a sequence.
- Melting Temperature Estimation: Estimates the melting temperature of a DNA sequence.
- FASTQ Filtering: Filters sequences based on GC content, length, and quality score.

Running from Terminal:
```bash
python main_dna_rna_filter.py
```

## Examples

1. Convert FASTA Multiline to Single Line::

```python
from bio_files_processor import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")
```

2. Parse BLAST Output::

```python
from bio_files_processor import parse_blast_output

parse_blast_output("blast_results.txt", "protein_list.txt")  
```

3. DNA/RNA Sequence Manipulation:

```python
from main_dna_rna_filter import DNASequence

dna_seq = DNASequence("ATGCGTACG")
print(dna_seq.transcribe().sequence) 
print(dna_seq.reverse_complement().sequence) 
```

4. FASTQ Filtering:

```python
from main_dna_rna_filter import FastqFilter

filter = FastqFilter(gc_bounds=(40, 60), length_bounds=(50, 200), quality_threshold=30)
filter.filter_fastq("input.fastq", "filtered.fastq")
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
Q: What types of sequences are supported?A: The scripts support DNA, RNA, and protein sequences.

Q: Can I filter sequences by GC content and quality score?A: Yes, you can set GC content and quality score thresholds using FastqFilter.

Q: How do I extract protein names from BLAST output?A: Use the parse_blast_output() function from bio_files_processor.py.


## Citation
If you use bioinformatics-utils in your research, please cite this repository and the following publication:

Lebedeva ES, (2024). bioinformatics-utils: A comprehensive toolkit for DNA and RNA sequence analysis.

## Contact

For any issues or questions, please report directly to the GitHub issue tracker, or contact the author via email at EkaterinaLebedeva2612@gmail.com.

# bioinformatics-utils

**bioinformatics-utils** is a  toolkit that provides essential utilities for working with biological sequences (DNA, RNA, proteins) and FASTQ files. It is designed for sequence filtering, manipulation, quality control, and parsing of BLAST results.

The repository includes:

- `bio_files_processor.py` — utility functions for processing FASTA files and parsing BLAST output  
- `main_dna_rna_filter.py` — core classes for biological sequence representation and FASTQ filtering (based on GC content, length, and quality)  
- `main.py` — command-line interface (CLI) for running FASTQ filtering with customizable arguments  
- `fastq_filter_test.py` — unit tests implemented using pytest to validate filtering logic, error handling, and file Input/Output  
- `logs/filter.log` — automatically generated log file storing filter warnings and errors  
- `examples/` — folder containing example input and output files, including sample `.fastq`, `.fasta`, `.txt` and BLAST files


## Authors:
- **Software Development**: Ekaterina Lebedeva, Bioinformatics Institute
- **Testing**: Ekaterina Lebedeva, Bioinformatics Institute
- **Supervisor**: Nikita Vaulin, Bioinformatics Institute

---


## Content
- [Installation](#installation)
- [Running Instructions](#running-instructions)
  - [FASTA Processing and BLAST Output Parsing](#fasta-processing-and-blast-output-parsing)
    - [Running from Terminal](#running-from-terminal)
    - [Running from Jupyter Notebooks or interactive Python environments](#running-from-jupyter-notebooks-or-interactive-python-environments)
  - [Sequence Manipulation & FASTQ Filtering (main_dna_rna_filter.py)](#sequence-manipulation--fastq-filtering-maindna_rna_filterpy)
    - [Running from Terminal](#running-from-terminal-1)
    - [Running from Jupyter Notebooks or interactive Python environments](#running-from-jupyter-notebooks-or-interactive-python-environments-1)
- [Running Tests](#running-tests)
- [FAQ](#faq)
- [Citation](#citation)
- [Contacts](#contact)
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


## Running Instructions
**bioinformatics-utils** offers two main functionalities:

### **FASTA Processing and BLAST Output Parsing:**
Functions:

- convert_multiline_fasta_to_oneline: Converts multi-line FASTA files into single-line format.
- parse_blast_output: Extracts the top protein name from each query in BLAST output files.

#### Running from Terminal:
```bash
python bio_files_processor.py
```
You will be prompted to choose a function and optionally specify input/output file paths.
If no input is given, default files from the examples/ folder will be used.
Follow the prompt to select:

- Enter 1 for FASTA conversion.
    - Enter path to input and output  FASTA files
```bash
Enter path to input FASTA file (default: examples/example_multiline_fasta.fasta): examples/example_multiline_fasta.fasta
Enter path to output FASTA file (default: examples/output_oneline_fasta.fasta): examples/output_oneline_fasta.fasta
```

- Enter 2 for BLAST output parsing
    - Enter path to input and output txt files
```bash
Enter path to BLAST output file (default: examples/example_blast_results.txt): examples/example_blast_results.txt
Enter path to output protein list (default: examples/protein_list.txt): examples/protein_list.txt
```
#### Running from Jupyter Notebooks or interactive Python environments

1. Convert FASTA Multiline to Single Line::

```python
from bio_files_processor import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")
```

2. Parse BLAST Output:

```python
from bio_files_processor import parse_blast_output

parse_blast_output("blast_results.txt", "protein_list.txt")  
```

### Sequence Manipulation & FASTQ Filtering (main_dna_rna_filter.py)

Features:

- Transcription: DNA → RNA conversion.
- Reverse Complement: Computes the reverse complement of a sequence.
- GC Content Calculation: Computes the GC content of a sequence.
- Melting Temperature Estimation: Estimates the melting temperature of a DNA sequence.
- FASTQ Filtering: Filters sequences based on GC content, length, and quality score.

All core functionalities of ``main_dna_rna_filter.py`` can be used in Jupyter Notebooks or interactive Python environments. 

### Running from Terminal
In addition, FASTQ filtering can also be run directly from the terminal using main.py with command-line arguments.
```bash
python main.py \
  -i examples/example_fastq.fastq \
  -o examples/output_filtered.fastq \
  --min_gc 40 \
  --max_gc 60 \
  --min_len 50 \
  --max_len 200 \
  --min_quality 30
```

| Argument           | Description                          |   Example                     |
|--------------------|--------------------------------------|-------------------------------|
| `-i`, `--input`    | Path to the input FASTQ file          | `-i examples/input.fastq`   |
| `-o`, `--output`   | Path to the output FASTQ file           | `-o examples/output.fastq`  |
| `--min_gc`         | Minimum GC content percentage         | `--min_gc 40`               |
| `--max_gc`         | Maximum GC content percentage          | `--max_gc 60`               |
| `--min_len`        | Minimum read length                    | `--min_len 50`              |
| `--max_len`        | Maximum read length                   | `--max_len 200`             |
| `--min_quality`    | Minimum average Phred quality score  | `--min_quality 30`          |

**Note:** All optional parameters have default values:
- `min_gc=0`, `max_gc=100`
- `min_len=0`, `max_len=4294967296`
- `min_quality=0`

#### Running from Jupyter Notebooks or interactive Python environments

1. DNA/RNA Sequence Manipulation:

```python
from main_dna_rna_filter import DNASequence

# Create a DNA sequence
dna_seq = DNASequence("ATGCGTACG")

# Transcribe DNA → RNA
rna_seq = dna_seq.transcribe()
print("RNA:", rna_seq.sequence)

# Get reverse complement
print("Reverse Complement:", dna_seq.reverse_complement().sequence)

# Calculate GC content
print("GC Content (%):", round(dna_seq.gc_content(), 2)) 

# Estimate melting temperature
print("Melting Temperature (°C):", dna_seq.melting_temperature()) 

# Use RNASequence methods
print("RNA Complement:", rna_seq.complement().sequence)
print("RNA Reverse:", rna_seq.reverse().sequence)

# Create a protein sequence and calculate molecular weight
protein = AminoAcidSequence("MVHLTPEEKSAVTAL")
print("Molecular Weight (Da):", round(protein.molecular_weight(), 2))
```
2. FASTQ Filtering:

```python
from main_dna_rna_filter import FastqFilter

# Initialize filtering tool with desired thresholds
filter = FastqFilter(
    gc_bounds=(40, 60), 
    length_bounds=(50, 200),   

# Run filtering on an input FASTQ file
filter.filter_fastq("input.fastq", "filtered.fastq")
```
## Running Tests
Unit tests are implemented using pytest.They validate core functionalities such as:

- FASTQ filtering logic (GC content, quality, read length)
- Error handling (e.g., missing quality scores, empty reads)
- File I/O (reading and writing FASTQ files)
- GC/quality edge cases

To run the tests, simply execute:

```bash
pytest fastq_filter_test.py
```
All tests are self-contained and cover both successful filtering and error scenarios. Temporary files are automatically created and cleaned up using tmp_path.

## FAQ
**Q: What types of sequences are supported?**  
**A:** The toolkit supports DNA, RNA, and protein sequences. You can manipulate them using dedicated classes like `DNASequence`, `RNASequence`, and `AminoAcidSequence`.

---

**Q: Can I filter sequences by GC content and quality score?**  
**A:** Yes. You can set GC content range, minimum read length, and average Phred quality score using the `FastqFilter` class or the CLI (`main.py`). Filtered results are saved to a new FASTQ file, and discarded records are logged in `logs/filter.log`.

---

**Q: How do I extract protein names from BLAST output?**  
**A:** Use the `parse_blast_output()` function from `bio_files_processor.py`. It parses tabular BLAST output and extracts the top protein name for each query.

---

**Q: Can I use this toolkit in Jupyter Notebooks?**  
**A:** Yes! All classes and functions can be imported and used interactively. Only `main.py` is designed for command-line usage.

---

**Q: Are example input files included?**  
**A:** Yes. You can find sample FASTQ, FASTA, and BLAST files in the `examples/` directory to test all major functions.


## Citation
If you use bioinformatics-utils in your research, please cite this repository and the following publication:

Lebedeva ES, (2024). bioinformatics-utils: A comprehensive toolkit for DNA and RNA sequence analysis.

## Contact

For any issues or questions, please report directly to the GitHub issue tracker, or contact the author via email at EkaterinaLebedeva2612@gmail.com.

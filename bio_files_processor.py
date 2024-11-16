from typing import Optional


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: Optional[str] = None
) -> None:
    """
    Converts a multi-line FASTA file to a single-line FASTA file where each
    sequence is on a single line.

    Arguments:
    - input_fasta: Path to the input multi-line FASTA file.
    - output_fasta: Path to the output single-line FASTA file. If not
      provided, the output is printed to the terminal.
    """
    output_file = open(output_fasta, 'w') if output_fasta else None

    with open(input_fasta, 'r') as file:
        current_header = None
        current_sequence = []

        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    sequence = ''.join(current_sequence)
                    if output_file:
                        output_file.write(current_header + "\n")
                        output_file.write(sequence + "\n")
                    else:
                        print(current_header)
                        print(sequence)
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_header:
            sequence = ''.join(current_sequence)
            if output_file:
                output_file.write(current_header + "\n")
                output_file.write(sequence + "\n")
            else:
                print(current_header)
                print(sequence)


def parse_blast_output(
    input_file: str, output_file: Optional[str] = None
) -> None:
    """
    Parses a BLAST output file and extracts the top protein name for each
    query. The results are sorted alphabetically and saved to the output
    file.

    Arguments:
    - input_file: Path to the BLAST output file (txt).
    - output_file: Path to the output file (optional). If not provided,
      results are printed.
    """
    results = []
    capture = False
    capture_description = False

    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith("Query #"):
                capture = False
                capture_description = False

            elif line.startswith(
                "Sequences producing significant alignments:"
            ):
                capture = True

            elif capture and line.startswith("Description"):
                capture_description = True
                continue

            elif capture_description and not line.startswith(" "):
                description = line.split("    ")[0]
                results.append(description)
                capture_description = False

    results.sort()

    if output_file:
        with open(output_file, 'w') as out_file:
            for result in results:
                out_file.write(result + "\n")
    else:
        for result in results:
            print(result)


if __name__ == "__main__":
    input_fasta = "example_multiline_fasta.fasta"
    output_fasta = "output.fasta"
    input_file = "example_blast_results.txt"
    output_file = "protein_list.txt"

    print("Выберите функцию для запуска:")
    print("1 - convert_multiline_fasta_to_oneline")
    print("2 - parse_blast_output")
    choice = input("Введите 1 или 2: ")

    if choice == "1":
        convert_multiline_fasta_to_oneline(input_fasta, output_fasta)
    elif choice == "2":
        parse_blast_output(input_file, output_file)
    else:
        print("Неверный выбор. Запустите скрипт снова.")

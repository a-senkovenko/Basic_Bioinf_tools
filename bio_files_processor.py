import os


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> None:
    """
    Convert multiline fasta to oneline
    Arguments:
        input_fasta (str): path to input FASTA
        output_fasta (str): name of output file
    Returns: 
        None. Write into new file with name 'output_fasta'
    """
    output_dir = "convert_results"
    os.makedirs(output_dir, exist_ok=True)

    if output_fasta is None:
        base_name = os.path.basename(input_fasta)
        name, ext = os.path.splitext(base_name)
        output_fasta = f"{name}_oneline{ext or '.fasta'}"

    output_path = os.path.join(output_dir, output_fasta)

    with open(input_fasta, "r") as inp, open(output_path, "w") as res:
        header = None
        seq_parts = []

        for line in inp:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    res.write(header + "\n")
                    res.write("".join(seq_parts) + "\n")
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)

        if header:
            res.write(header + "\n")
            res.write("".join(seq_parts) + "\n")


def parse_blast_output():
    pass

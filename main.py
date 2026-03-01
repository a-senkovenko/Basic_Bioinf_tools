import os
from typing import Iterator, Tuple, List, Dict

MAX_GC = 100.0
MAX_LENGTH = 2**32
DNA_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
RNA_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G", "N": "N"}
TRANS_MAP = {"T": "U", "t": "u", "U": "T", "u": "t"}

def read_fastq(path: str) -> Iterator[Tuple[str, str, str]]:
    """
    Read FASTQ file and make an iterator(header, seq, qual).
    """
    with open(path, "r") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()

            if not (header.startswith("@") and plus.startswith("+")):
                raise ValueError(f"Invalid FASTQ format near {header}")

            yield (header, seq, qual)


def write_fastq(sequences: Iterator[Tuple[str, str, str]],
                output_path: str) -> None:
    """
    Write iterator(header, seq, qual) to file.
    Make folder if it doesn't exist.
    Overwrite is disabled.
    """
    output_path = os.path.join("filtered", output_path)
    folder = os.path.dirname(output_path)
    if folder and not os.path.exists(folder):
        os.makedirs(folder)

    if os.path.exists(output_path):
        raise FileExistsError(f"File {output_path} already exists.")

    with open(output_path, "w") as f:
        for header, seq, qual in sequences:
            f.write(f"@{header}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{qual}\n")


def compute_gc_content(sequence: str) -> float:
    """Return GC percentage of given sequence."""
    if not sequence:
        return 0.0
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


def quality_score(quality: str) -> float:
    """Return average quality score (Phred+33) of given sequence."""
    if not quality:
        return 0.0
    return sum(ord(ch) - 33 for ch in quality) / len(quality)


def length_filter(
    seqs: Dict[str, Tuple[str, str]],
    length_bounds: Tuple[int, int] = (0, MAX_LENGTH)
) -> List[str]:
    """Return names of reads filtered by length."""
    min_len, max_len = length_bounds
    return [name for name, (seq, _) in seqs.items()
            if min_len <= len(seq) <= max_len]


def gc_filter(
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Tuple[float, float] = (0.0, MAX_GC)
) -> List[str]:
    """Return names of reads filtered by GC percentage."""
    min_gc, max_gc = gc_bounds
    return [
        name
        for name, (seq, _) in seqs.items()
        if min_gc <= compute_gc_content(seq) <= max_gc
    ]


def quality_filter(
    seqs: Dict[str, Tuple[str, str]], quality_threshold: float = 0.0
) -> List[str]:
    """Return names of reads with average quality above threshold."""
    return [
        name
        for name, (_, qual) in seqs.items()
        if quality_score(qual) >= quality_threshold
    ]


def filter_fastq(
    seqs: Iterator[Tuple[str, str, str]],
    gc_bounds: Tuple[float, float] | float = (0.0, MAX_GC),
    length_bounds: Tuple[int, int] | int = (0, MAX_LENGTH),
    quality_threshold: float = 0.0,
) -> Iterator[Tuple[str, str, str]]:
    """
    Filter FASTQ sequences by GC content, length, and quality.
    Uses iteration for fast filtration by lines in file.
    Arguments:
        seqs (Dict[str, Tuple[str, str]]): Dictionary of sequences
        gc_bounds (Tuple[float, float]): GC composition (%)
        length_bounds (Tuple[int, int]): length for filtration
        quality_threshold (float): Phred33 scale quality, default = 0
    Returns:
        Iterator[Tuple[str, str, str]] for function run_filter_fastq()
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0.0, float(gc_bounds))

    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    for header, seq, qual in seqs:
        if (
            length_filter({header: (seq, qual)}, length_bounds)
            and gc_filter({header: (seq, qual)}, gc_bounds)
            and quality_filter({header: (seq, qual)}, quality_threshold)
        ):
            yield (header, seq, qual)


def run_dna_rna_tools(*args: str) -> List[str]:
    """
    Check input format. Args must be 'str' type, positional.
    The last arg must be procedure.
    If wrong format: ValueError.
    Run procedures and collect results as list.
    If unknown procedure: KeyError.
    Arguments:
        seq: Sequences to run procedures.
        procedure: operation to perform on seq.
    Returns:
        Results of run procedure as list.
    """

    if len(args) < 2:
        raise ValueError("Need at least one sequence and procedure")
    *sequences, procedure = args

    procedures = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    results = []
    for seq in sequences:
        if procedure not in procedures:
            raise KeyError(f"Unknown procedure: {tuple(procedures.keys())}")
        if procedure != "is_nucleic_acid":
            if not is_nucleic_acid(seq):
                raise ValueError(f"Wrong or ambiguous sequence: {seq}")
        result = procedures[procedure](seq)
        results.append(result)

    if len(results) == 1:
        return results[0]
    return results


def run_filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0.0, MAX_GC),
    length_bounds=(0, MAX_LENGTH),
    quality_threshold=0.0,
) -> None:
    """
    Fast filter FASTQ sequences by GC content, length, and quality.
    Uses iteration for fast filtration by lines in file.
    Arguments:
        input_fastq (str): name of input file
        output_fastq (str): name of output file
        gc_bounds (Tuple[float, float]): GC composition (%)
        length_bounds (Tuple[int, int]): length for filtration
        quality_threshold (float): Phred33 scale quality, default = 0
    Returns:
        None. Write filtered sequences to file with name 'output_fastq'
    """
    filtered = filter_fastq(
        read_fastq(input_fastq),
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=quality_threshold,
    )
    write_fastq(filtered, output_fastq)
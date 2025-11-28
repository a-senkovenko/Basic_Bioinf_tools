from modules.dna_rna_tools import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement,
)
from modules.filter_fastq import filter_fastq, MAX_GC, MAX_LENGTH
from modules.io_fastq import read_fastq, write_fastq
from typing import List


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

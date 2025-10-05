import sys

import modules.dna_rna_tools as drt
import modules.filter_fastq as ff
from typing import Dict, Tuple, Union


def run_dna_rna_tools(*args):
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
        result = procedures[procedure](seq)
        results.append(result)

    if len(results) == 1:
        return results[0]
    return results


def filter_fastq(
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Tuple[float, float] = (0.0, MAX_GC),
    length_bounds: Tuple[int, int] = (0, MAX_LENGTH),
    quality_threshold: float = 0.0,
) -> Dict[str, Tuple[str, str]]:
    """
    Filter FASTQ sequences by GC content, length, and quality.
    Check type of 'seqs' argument.
    Wrong 'seqs': TypeError.
    Arguments:
        seqs (Dict[str, Tuple[str, str]]): Dictionary of sequences
        gc_bounds (Tuple[float, float]): GC composition (%)
        length_bounds (Tuple[int, int]): length for filtration
        quality_threshold (float): Phred33 scale quality, default = 0
    Returns:
        result (Dict[str, Tuple[str, str]]): Filtered sequences
    """
    length_pass = set(length_filter(seqs, length_bounds))
    gc_pass = set(gc_filter(seqs, gc_bounds))
    quality_pass = set(quality_filter(seqs, quality_threshold))

    valid_names = length_pass & gc_pass & quality_pass

    result = {name: seqs[name] for name in valid_names}

    return result

DNA_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
RNA_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G", "N": "N"}
TRANS_MAP = {"T": "U", "t": "u", "U": "T", "u": "t"}


def is_nucleic_acid(seq: str) -> bool:
    """
    Check if sequence is nucleic acid.
    Args:
        seq (str): Sequence to check.
    Returns:
        bool: True if sequence is DNA or RNA, otherwise False.
    """
    bases = set(seq.upper())
    valid_chars_dna = {"A", "G", "C", "T"}
    valid_chars_rna = {"A", "G", "C", "U"}

    return bases <= valid_chars_dna or bases <= valid_chars_rna


def transcribe(seq: str) -> str:
    """
    Apply transcription operation for nucleic acid sequences
    (reverse transcription for RNA sequences),
    check if nucleic acid sequences before operation with is_nucleic_acid().
    If not nucleic acid: ValueError.
    If seq has T and U: ValueError.
    Args:
        seq (str): Sequence to transcribe.
    Returns:
        seq (str): "ATGC" -> "AUGC"; "AUGC" -> "ATGC".
    """
    if not is_nucleic_acid(seq):
        raise ValueError(f"Wrong sequence: {seq}")
    if "T" in seq.upper() and "U" in seq.upper():
        raise ValueError(f"Ambiguous sequence: contains both T and U: {seq}")

    return "".join(TRANS_MAP.get(nucl, nucl) for nucl in seq)


def reverse(seq: str) -> str:
    """
    Reverts sequence, check sequences before operation with is_nucleic_acid().
    If not nucleic acid: ValueError.
    Args:
        seq (str): Sequence to revert.
    Returns:
        seq (str): "ATGC" -> "CGTA".
    """
    if not is_nucleic_acid(seq):
        raise ValueError(f"Wrong sequence: {seq}")

    return seq[::-1]


def complement(seq: str) -> str:
    """
    Build complementary, check before operation with is_nucleic_acid().
    Can build complementary for RNA sequences.
    If not nucleic acid: ValueError.
    Args:
        seq (str): Matrix sequence for complementary sequence.
    Returns:
        seq (str): "ATGC" -> "TACG"; "AUGC" -> "UACG".
    """
    if not is_nucleic_acid(seq):
        raise ValueError(f"Wrong sequence: {seq}")
    if "T" in seq.upper() and "U" in seq.upper():
        raise ValueError(f"Ambiguous sequence: contains both T and U: {seq}")

    complement_map = DNA_COMPLEMENT if "T" in seq.upper() else RNA_COMPLEMENT
    return "".join(complement_map.get(nucl, nucl) for nucl in seq)


def reverse_complement(seq: str) -> str:
    """
    Build complementary sequence and revert it,
    use complement() and reverse() functions consecutively.
    Can build complementary for RNA sequences.
    If not nucleic acid: ValueError.
    Args:
        seq: Matrix sequence for reverted complementary sequence.
    Returns:
        Complement sequence as a string, f.e. input: "ATGC" -> output: "GCAT".
    """
    return reverse(complement(seq))
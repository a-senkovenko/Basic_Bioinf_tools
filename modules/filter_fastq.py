from typing import Dict, Tuple, List

MAX_GC = 100.0
MAX_LENGTH = 2**32


def compute_gc_content(sequence: str) -> float:
    """Return GC percentage of given sequence."""
    if not sequence:
        return 0.0
    sequence = sequence.upper()
    gc_count = sequence.count("G")+sequence.count("C")
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
    return [
        name
        for name, (seq, _) in seqs.items()
        if min_len <= len(seq) <= max_len
        ]


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

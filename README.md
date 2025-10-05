# Basic Bioinf tools

This repository contains simple tools to perform basic operations with DNA/RNA sequences and FASTQ-format data.

## Content

The script was made with Python 3.11.

Two modules are available now:

1.  **Toolbox for DNA and RNA operations (dna_rna_tools.py)**:
    *   Input arguments must be strings. Needs at least one sequence to analyze and one procedure to run as last argument.
    *   Performs one of the following procedures:
        *   `is_nucleic_acid`: Check if sequence is nucleic acid.
        *   `transcribe`: Perform transcription and reverse transcription of given sequence.
        *   `reverse`: Reverse the sequence.
        *   `complement`: Build complementary sequence.
        *   `reverse_complement`: Build reversed complementary sequence.
    *   Returns result of given procedure as list of strings, saves register of input sequences.

2.  **Tool for FASTQ format filtering (filter_fastq.py)**:
    *   Input arguments must be in following format: 
        * seqs: Dict[str, Tuple[str, str]] - sequences to filter
        * gc_bounds: Tuple[float, float] - bounds for filtering by GC percentage
        * length_bounds: Tuple[int, int] - bounds for filtering by length
        * quality_threshold: float - threshold for filtering by quality
    *   Reads sequences in FASTQ format and filters the by given parameters.
    *   Implemented functions:
        *   `compute_gc_content`: Return GC percentage of given sequence.
        *   `quality_score`: Return average quality score (Phred+33) of given sequence.
        *   `length_filter`: Return names of reads filtered by length.
        *   `gc_filter`:Return names of reads filtered by GC percentage.
        *   `quality_filter`:Return names of reads with average quality above threshold.
    *   Returns filtered dictionary of the same structure as argument 'seqs'.

## Adjustments
    *   Add global constants in modules
    *   Rebuild `is_nucleic_acid` by reviewer's idea
    *   Fix `transcribe`

## TODO
    *   Add error handling in filter_fastq.py
    *   Add operations with files for filter_fastq.py


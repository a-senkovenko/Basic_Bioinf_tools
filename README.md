# Basic Bioinf tools

This repository contains simple tools to perform basic operations with DNA/RNA sequences and FASTQ-format data.

## Content

The script was made with Python 3.11.

Three modules are available now:

1.  **Toolbox for DNA and RNA operations (dna_rna_tools.py)**:
    *   Input arguments must be strings. Needs at least one sequence to analyze and one procedure to run as last argument.
    *   Performs one of the following procedures:
        *   `is_nucleic_acid`: Check if sequence is nucleic acid.
        *   `transcribe`: Perform transcription and reverse transcription of given sequence.
        *   `reverse`: Reverse the sequence.
        *   `complement`: Build complementary sequence.
        *   `reverse_complement`: Build reversed complementary sequence.
    *   Returns result of given procedure as list of strings, saves register of input sequences.

2.  **Improved tool for FASTQ format filtering (filter_fastq.py)**:
    *   Input arguments must be in following format: 
        * input_fastq (str): name of input file
        * output_fastq (str): name of output file
        * gc_bounds: Tuple[float, float] - bounds for filtering by GC percentage
        * length_bounds: Tuple[int, int] - bounds for filtering by length
        * quality_threshold: float - threshold for filtering by quality
    *   Reads sequences in FASTQ format from file and filters by the given parameters.
    *   Implemented functions:
        *   `compute_gc_content`: Return GC percentage of given sequence.
        *   `quality_score`: Return average quality score (Phred+33) of given sequence.
        *   `length_filter`: Return names of reads filtered by length.
        *   `gc_filter`:Return names of reads filtered by GC percentage.
        *   `quality_filter`:Return names of reads with average quality above threshold.
        *   `filter_fastq`:Return iterator(name, seq, qual).
    *   Writes down the file 'output_fastq' in folder 'filtered' with overwrite protection.
3.  **New module io_fastq.py**:
    *   Implemented functions:
        * `read_fastq`: Read FASTQ file and make an iterator(name, seq, qual) for 'filter_fastq()'
        * `write_fastq`: Write iterator(name, seq, qual) to file after filtration.


## Updates (05.10.2025)
    *   Add global constants in modules
    *   Rebuild `is_nucleic_acid` by reviewer's idea
    *   Fix `transcribe`

## Updates (12.10.2025)
    *   Fix data types in function definitions
    *   Fix repeating code in dna_rna_tools.py
    *   Move filter_fastq() from main.py into filter_fastq.py
    *   New run_filter_fastq() in main.py to operate filter_fastq.py quickly
    *   Add module io_fastq to control file input and output for run_filter_fastq()
    *   Small fixes in filter_fastq.py
    *   Fix imports in main.py
    *   Add new tool to work with FASTA format and BLAST output -> bio_files_processor.py
    *   Add new function -> convert_multiline_fasta_to_oneline() in bio_files_processor.py

## TODO
    *   Add parse_blast_output() into bio_files_processor.py


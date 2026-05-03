import os
import click
import logging, logger
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from typing import Tuple


logger = logging.getLogger(__name__)

class BiologicalSequence:
    """
    Attributes:
        seq: Sequences to run procedures.
        procedure: operation to perform on seq.
    """
    def __init__(self, seq:str):
        self._sequence = seq.upper()

    def __len__(self):
        return len(self._sequence)
    
    def __getitem__(self, item):
        return self._sequence[item]
    
    def _check_alphabet(self):
        """
        Compare with the given class alphabet.
        """
        if not set(self._sequence) <= self._alphabet:
            raise ValueError("Invalid nucleic acid sequence")

    def __repr__(self):
        return f"{self.__class__.__name__}('{self._sequence}')"

class NucleicAcidSequence(BiologicalSequence):
    """
    Class to store and operate nucleic acid sequence.
    Attributes:
        seq(str): Sequence to run procedures.
    """
    _complement = {}
    _alphabet = set()

    def complement(self):
        """
        Build complementary, check before operation with is_nucleic_acid().
        Can build complementary for RNA sequences.
        If not nucleic acid: ValueError.
        Args:
            seq (str): Matrix sequence for complementary sequence.
        Returns:
            seq (str): "ATGC" -> "TACG"; "AUGC" -> "UACG".
        """
        return "".join(self._complement[n] for n in self._sequence)

    def reverse(self):
        """
        Reverts sequence, check sequences before operation with is_nucleic_acid().
        If not nucleic acid: ValueError.
        Args:
            seq (str): Sequence to revert.
        Returns:
            seq (str): "ATGC" -> "CGTA".
        """
        return self[::-1]

    def reverse_complement(self):
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
        return self.complement()[::-1]

class DNASequence(NucleicAcidSequence):
    """
    """
    _alphabet = {"A", "T", "G", "C"}
    _complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    _trans_map = {"A": "A", "T": "U", "G": "G", "C": "C"}

    def transcribe(self):
        """
        Apply transcription operation for DNA.
        Args:
            seq (str): Sequence to transcribe.
        Returns:
            seq (str): "ATGC" -> "AUGC".
        """
        return "".join(self._trans_map[n] for n in self._sequence)

class RNASequence(NucleicAcidSequence):
    _alphabet = {"A", "U", "G", "C"}
    _complement = {"A": "U", "U": "A", "G": "C", "C": "G"}
    _trans_map = {"A": "A", "U": "T", "G": "G", "C": "C"}

    def reverse_transcribe(self):
        """
        Apply transcription operation for RNA.
        Args:
            seq (str): Sequence to transcribe.
        Returns:
            seq (str): "AUGC" -> "ATGC".
        """
        return "".join(self._trans_map[n] for n in self._sequence)

class AminoAcidSequence(BiologicalSequence):
    _alphabet = {
        "A", "R", "N", "D", "C", "Q", "E", "G",
        "H", "I", "L", "K", "M", "F", "P", "S",
        "T", "W", "Y", "V"
    }

    def molecular_weight(self):
        weights = {
            "A": 89, "R": 174, "N": 132, "D": 133, "C": 121,
            "Q": 146, "E": 147, "G": 75, "H": 155, "I": 131,
            "L": 131, "K": 146, "M": 149, "F": 165, "P": 115,
            "S": 105, "T": 119, "W": 204, "Y": 181, "V": 117
        }
        return sum(weights[aa] for aa in self._sequence)

def setup_logging():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    )

    info_handler = logging.FileHandler("info.log")
    info_handler.setLevel(logging.INFO)

    class InfoFilter(logging.Filter):
        def filter(self, record):
            return record.levelno < logging.ERROR

    info_handler.addFilter(InfoFilter())
    info_handler.setFormatter(formatter)

    error_handler = logging.FileHandler("error.log")
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(formatter)

    logger.addHandler(info_handler)
    logger.addHandler(error_handler)

def filter_fastq(
    input_path: str,
    output_path: str,
    gc_bounds: Tuple[float, float] | float = (0.0, 100.0),
    length_bounds: Tuple[int, int] | int = (0, float("inf")),
    quality_threshold: float = 0.0
) -> None:
    """
    Filter FASTQ function with Biopython.

    Args:
        input_path: input path for FASTQ
        output_path: output path for FASTQ
        gc_bounds: GC% bounds
        length_bounds: length bounds
        quality_threshold: mean Phred quality
    """
    logger.info("Starting FASTQ filtering")
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_path}")

    try:
        if isinstance(gc_bounds, (int, float)):
            gc_bounds = (0.0, float(gc_bounds))

        if isinstance(length_bounds, int):
            length_bounds = (0, length_bounds)

        min_gc, max_gc = gc_bounds
        min_len, max_len = length_bounds

        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

        if os.path.exists(output_path):
            logger.error("Output file already exists")
            raise FileExistsError(f"{output_path} already exists")

        filtered_records = []

        for record in SeqIO.parse(input_path, "fastq"):
            seq = record.seq
            length = len(seq)

            if not (min_len <= length <= max_len):
                continue

            gc = gc_fraction(seq) * 100
            if not (min_gc <= gc <= max_gc):
                continue

            qualities = record.letter_annotations["phred_quality"]
            avg_quality = sum(qualities) / len(qualities)

            if avg_quality < quality_threshold:
                continue

            filtered_records.append(record)

        SeqIO.write(filtered_records, output_path, "fastq")
    
        logger.info(f"Filtering finished successfully")

    except:
        logger.exception("Error during FASTQ filtering")
        raise

###CLI with click###

def parse_bounds(value: str, is_float=True) -> Tuple:
    if "," in value:
        a, b = value.split(",")
        return (float(a), float(b)) if is_float else (int(a), int(b))
    else:
        return float(value) if is_float else int(value)


@click.command()
@click.argument("input")
@click.argument("output")
@click.option("--gc", default="0,100", help="GC bounds (e.g. 30,70 or 50)")
@click.option("--length", default="0,inf", help="Length bounds")
@click.option("--quality", default=0.0, type=float, help="Min mean quality")
def cli(input, output, gc, length, quality):
    """Filter FASTQ file"""

    gc_bounds = parse_bounds(gc, is_float=True)

    if length == "inf":
        length_bounds = (0, float("inf"))
    else:
        length_bounds = parse_bounds(length, is_float=False)

    filter_fastq(
        input_path=input,
        output_path=output,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=quality
    )


if __name__ == "__main__":
    cli()
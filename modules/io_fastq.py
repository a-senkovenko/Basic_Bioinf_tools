import os
from typing import Iterator, Tuple


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

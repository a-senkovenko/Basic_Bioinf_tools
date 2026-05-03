import os
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from main import filter_fastq


def create_fastq(records, path):
    SeqIO.write(records, path, "fastq")


@pytest.fixture
def sample_records():
    """Generate test sequence data"""
    rec1 = SeqRecord(Seq("ATGC"), id="r1")
    rec1.letter_annotations["phred_quality"] = [40, 40, 40, 40]

    rec2 = SeqRecord(Seq("AAAA"), id="r2")
    rec2.letter_annotations["phred_quality"] = [10, 10, 10, 10]

    rec3 = SeqRecord(Seq("GGGGCCCC"), id="r3")
    rec3.letter_annotations["phred_quality"] = [30] * 8

    return [rec1, rec2, rec3]


def test_basic_filtering(tmp_path, sample_records):
    """Filtering test to check how the function work"""
    input_file = tmp_path / "in.fastq"
    output_file = tmp_path / "out.fastq"

    create_fastq(sample_records, input_file)

    filter_fastq(str(input_file), str(output_file), quality_threshold=20)

    result = list(SeqIO.parse(output_file, "fastq"))
    assert len(result) == 2
    assert {r.id for r in result} == {"r1", "r3"}


def test_gc_filter(tmp_path, sample_records):
    """GC-bound filtering"""
    input_file = tmp_path / "in.fastq"
    output_file = tmp_path / "out.fastq"

    create_fastq(sample_records, input_file)

    filter_fastq(str(input_file), str(output_file), gc_bounds=(50, 100))

    result = list(SeqIO.parse(output_file, "fastq"))
    assert all(50 <= (len([b for b in r.seq if b in "GC"]) / len(r.seq) * 100) <= 100 for r in result)


def test_length_filter(tmp_path, sample_records):
    """Len-bound filtering"""
    input_file = tmp_path / "in.fastq"
    output_file = tmp_path / "out.fastq"

    create_fastq(sample_records, input_file)

    filter_fastq(str(input_file), str(output_file), length_bounds=(5, 10))

    result = list(SeqIO.parse(output_file, "fastq"))
    assert all(5 <= len(r.seq) <= 10 for r in result)


def test_quality_filter(tmp_path, sample_records):
    """QC-bound filtering"""
    input_file = tmp_path / "in.fastq"
    output_file = tmp_path / "out.fastq"

    create_fastq(sample_records, input_file)

    filter_fastq(str(input_file), str(output_file), quality_threshold=35)

    result = list(SeqIO.parse(output_file, "fastq"))
    assert len(result) == 1
    assert result[0].id == "r1"


def test_output_file_created(tmp_path, sample_records):
    """Check output file"""
    input_file = tmp_path / "in.fastq"
    output_file = tmp_path / "out.fastq"

    create_fastq(sample_records, input_file)

    filter_fastq(str(input_file), str(output_file))

    assert os.path.exists(output_file)


def test_output_file_not_overwritten(tmp_path, sample_records):
    """check output file existence"""
    input_file = tmp_path / "in.fastq"
    output_file = tmp_path / "out.fastq"

    create_fastq(sample_records, input_file)
    create_fastq(sample_records, output_file)

    with pytest.raises(FileExistsError):
        filter_fastq(str(input_file), str(output_file))


def test_empty_result(tmp_path, sample_records):
    """Check empty output"""
    input_file = tmp_path / "in.fastq"
    output_file = tmp_path / "out.fastq"

    create_fastq(sample_records, input_file)

    filter_fastq(str(input_file), str(output_file), quality_threshold=100)

    result = list(SeqIO.parse(output_file, "fastq"))
    assert result == []


def test_single_value_bounds(tmp_path, sample_records):
    """Проверка передачи gc_bounds и length_bounds как одного числа"""
    input_file = tmp_path / "in.fastq"
    output_file = tmp_path / "out.fastq"

    create_fastq(sample_records, input_file)

    filter_fastq(
        str(input_file),
        str(output_file),
        gc_bounds=50,
        length_bounds=4
    )

    result = list(SeqIO.parse(output_file, "fastq"))

    for r in result:
        assert len(r.seq) <= 4
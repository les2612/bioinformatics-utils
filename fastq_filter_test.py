import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from main_dna_rna_filter import FastqFilter


@pytest.fixture
def sample_records():
    return [
        SeqRecord(
            Seq("ATGC"),
            id="good1",
            letter_annotations={"phred_quality": [40, 40, 40, 40]}
            ),
        SeqRecord(
            Seq("AAAA"),
            id="low_gc",
            letter_annotations={"phred_quality": [30, 30, 30, 30]}
            ),
        SeqRecord(
            Seq("ATGC"),
            id="low_qual",
            letter_annotations={"phred_quality": [5, 5, 5, 5]}
            ),
        SeqRecord(
            Seq(""),
            id="empty",
            letter_annotations={"phred_quality": []}
            ),
    ]


@pytest.fixture
def write_fastq(tmp_path):
    def _write(records, filename="input.fastq"):
        path = tmp_path / filename
        SeqIO.write(records, path, "fastq")
        return path
    return _write


class TestFiltering:
    def test_quality_filter(self, sample_records):
        f = FastqFilter(quality_threshold=10)
        result = f.filter_records([sample_records[2]])
        assert len(result) == 0

    def test_gc_filter(self, sample_records):
        f = FastqFilter(gc_bounds=(40, 60))
        result = f.filter_records([sample_records[1]])
        assert len(result) == 0

    def test_passes_all_filters(self, sample_records):
        f = FastqFilter(gc_bounds=(40, 60), quality_threshold=10)
        result = f.filter_records([sample_records[0]])
        assert len(result) == 1
        assert result[0].id == "good1"

    def test_short_sequence(self):
        short = SeqRecord(
            Seq("A"),
            id="too_short",
            letter_annotations={"phred_quality": [30]}
            )
        f = FastqFilter(length_bounds=(2, 1000))
        result = f.filter_records([short])
        assert len(result) == 0


class TestIO:
    def test_read_write_fastq(self, sample_records, write_fastq, tmp_path):
        input_fh = write_fastq([sample_records[0], sample_records[1]])
        output_fh = tmp_path / "output.fastq"

        f = FastqFilter(gc_bounds=(40, 60), quality_threshold=10)
        f.filter_fastq(str(input_fh), str(output_fh))

        result = list(SeqIO.parse(str(output_fh), "fastq"))
        assert len(result) == 1
        assert result[0].id == "good1"

    def test_missing_phred_quality(self, capfd):
        record = SeqRecord(Seq("ATGC"), id="no_qual")
        f = FastqFilter()
        f.filter_records([record])
        out, _ = capfd.readouterr()
        assert "missing phred" in out.lower()

    def test_empty_input_fastq_raises(self, tmp_path):
        input_fh = tmp_path / "empty.fastq"
        input_fh.write_text("")
        output_fh = tmp_path / "output.fastq"
        f = FastqFilter()
        with pytest.raises(ValueError, match="No records found"):
            f.filter_fastq(str(input_fh), str(output_fh))


class TestUtils:
    @pytest.mark.parametrize("qualities, expected", [
        ([10, 20, 30], 20.0),
        ([0, 0, 0], 0.0),
        ([40], 40.0),
        ([5, 5, 5, 5], 5.0)
    ])
    def test_avg_quality_calculation_param(self, qualities, expected):
        assert FastqFilter.avg_quality(qualities) == expected


class TestErrors:
    def test_invalid_gc_bounds_raises(self):
        with pytest.raises(ValueError):
            FastqFilter(gc_bounds=(80, 20))

    def test_empty_sequence_prints(self, sample_records, capfd):
        f = FastqFilter()
        f.filter_records([sample_records[3]])
        out, _ = capfd.readouterr()
        assert "empty" in out.lower()

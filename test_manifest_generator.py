import re
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from manifest_generator import Fastq, ManifestParser

# ─────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────

FASTQ_EXT = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

R1_PATTERN = re.compile(r"(?:_R?1)(?:_[\w]+)?$")
R2_PATTERN = re.compile(r"(?:_R?2)(?:_[\w]+)?$")

TEST_DIR = Path("test_dir")

SHORT_R1 = TEST_DIR / "sampleA_R1.fastq.gz"
SHORT_R2 = TEST_DIR / "sampleA_R2.fastq.gz"
SHORT_R1_S2 = TEST_DIR / "sampleB_R1.fastq.gz"
SHORT_R2_S2 = TEST_DIR / "sampleB_R2.fastq.gz"

LONG_1 = TEST_DIR / "sampleA.fastq.gz"
LONG_2 = TEST_DIR / "sampleB.fastq.gz"

HYBRID_R1 = TEST_DIR / "sampleA_R1.fastq.gz"
HYBRID_R2 = TEST_DIR / "sampleA_R2.fastq.gz"
HYBRID_LONG = TEST_DIR / "sampleA_long.fastq.gz"


@pytest.fixture
def short_df():
    return pd.DataFrame(
        {
            "ID": ["sample_1", "sample_1", "sample_2", "sample_2"],
            "R1": [SHORT_R1, None, SHORT_R1_S2, None],
            "R2": [None, SHORT_R2, None, SHORT_R2_S2],
        }
    )


@pytest.fixture
def long_df():
    return pd.DataFrame(
        {
            "ID": ["sample_1", "sample_2"],
            "long_fastq": [LONG_1, LONG_2],
            "genome_size": [pd.NA, pd.NA],
        }
    )


@pytest.fixture
def hybrid_df():
    return pd.DataFrame(
        {
            "ID": ["sample_1", "sample_1", "sample_1"],
            "R1": [HYBRID_R1, None, None],
            "R2": [None, HYBRID_R2, None],
            "long_fastq": [None, None, HYBRID_LONG],
            "genome_size": [pd.NA, pd.NA, pd.NA],
        }
    )


# ─────────────────────────────────────────────
# Fastq class
# ─────────────────────────────────────────────


class TestFastq:
    def test_get_id_r1(self):
        f = Fastq(path=SHORT_R1)
        f._get_id()
        assert f.ID == "sampleA"

    def test_get_id_r2(self):
        f = Fastq(SHORT_R2)
        f._get_id()
        assert f.ID == "sampleA"

    def test_get_id_long(self):
        f = Fastq(LONG_1)
        f._get_id()
        assert f.ID == "sampleA"

    def test_detect_read_type_r1(self):
        assert Fastq(SHORT_R1)._detect_read_type() == "R1"

    def test_detect_read_type_r2(self):
        assert Fastq(SHORT_R2)._detect_read_type() == "R2"

    def test_detect_read_type_long(self):
        assert Fastq(LONG_1)._detect_read_type() == "long_fastq"

    def test_assign_properties_short_r1(self):
        f = Fastq(SHORT_R1_S2)
        f.assign_properties("short")
        assert f.R1 == SHORT_R1_S2
        assert f.R2 is None
        assert f.long_fastq is None
        assert f.ID == "sampleB"

    def test_assign_properties_long(self):
        f = Fastq(LONG_2)
        f.assign_properties("long")
        assert f.long_fastq == LONG_2
        assert f.R1 is None
        assert f.R2 is None
        assert f.ID == "sampleB"

    def test_assign_properties_hybrid_skips_id(self):
        """In hybrid mode, ID is pre-assigned and should not be overwritten"""
        f = Fastq(path=HYBRID_R1, ID="sample_1")
        f.assign_properties("hybrid")
        assert f.ID == "sample_1"
        assert f.R1 == HYBRID_R1

    def test_to_dict_keys(self):
        f = Fastq(SHORT_R1)
        d = f.to_dict()
        assert set(d.keys()) == {"ID", "R1", "R2", "long_fastq", "genome_size"}


# ─────────────────────────────────────────────
# ManifestParser — duplicate detection
# ─────────────────────────────────────────────


class TestManifestParserDuplicateDetection:
    def test_no_duplicates_long(self, long_df):
        parser = ManifestParser(long_df, "long", "error")
        assert parser._detect_duplicates() == (False, 0)

    def test_no_duplicates_short(self, short_df):
        parser = ManifestParser(short_df, "short", "error")
        assert parser._detect_duplicates() == (False, 0)

    def test_no_duplicates_hybrid(self, hybrid_df):
        parser = ManifestParser(hybrid_df, "hybrid", "error")
        assert parser._detect_duplicates() == (False, 0)

    def test_detects_duplicate_long(self):
        df = pd.DataFrame(
            {
                "ID": ["sample_1", "sample_1"],
                "long_fastq": [LONG_1, LONG_1],
                "genome_size": [pd.NA, pd.NA],
            }
        )
        parser = ManifestParser(df, "long", "warn")
        has_dup, n = parser._detect_duplicates()
        assert has_dup is True
        assert n > 0

    def test_detects_duplicate_short(self):
        df = pd.DataFrame(
            {
                "ID": ["sample_1", "sample_1", "sample_1"],
                "R1": [SHORT_R1, SHORT_R1, None],
                "R2": [None, None, SHORT_R2],
            }
        )
        parser = ManifestParser(df, "short", "warn")
        has_dup, n = parser._detect_duplicates()
        assert has_dup is True
        assert n > 0


# ─────────────────────────────────────────────
# ManifestParser — handle_duplicates methods
# ─────────────────────────────────────────────


class TestManifestParserHandleDuplicates:
    def test_error_exits_on_duplicate(self):
        df = pd.DataFrame(
            {
                "ID": ["sample_1", "sample_1", "sample_1"],
                "long_fastq": [LONG_1, LONG_1, LONG_1],
                "genome_size": [pd.NA, pd.NA, pd.NA],
            }
        )
        parser = ManifestParser(df=df, mode="long", method="error")
        with pytest.raises(SystemExit):
            parser.handle_duplicates()

    def test_warn_returns_df(self, short_df):
        dup_df = pd.concat([short_df, short_df])
        parser = ManifestParser(df=dup_df, mode="short", method="warn")
        result = parser.handle_duplicates()
        assert isinstance(result, pd.DataFrame)

    def test_drop_deduplicates(self, short_df):
        dup_df = pd.concat([short_df, short_df])
        parser = ManifestParser(df=dup_df, mode="short", method="drop")
        result = parser.handle_duplicates()
        assert result["ID"].duplicated().sum() == 0

    def test_no_duplicate_short_collapses(self, short_df):
        parser = ManifestParser(df=short_df, mode="short", method="error")
        result = parser.handle_duplicates()
        assert "sample_1" in result["ID"].values
        assert list(result.columns) == ["ID", "R1", "R2"]

    def test_no_duplicate_long_collapses(self, long_df):
        parser = ManifestParser(df=long_df, mode="long", method="error")
        result = parser.handle_duplicates()
        assert len(result) == 2

    def test_no_duplicate_hybrid_collapses(self, hybrid_df):
        parser = ManifestParser(df=hybrid_df, mode="hybrid", method="error")
        result = parser.handle_duplicates()
        assert "sample_1" in result["ID"].values


# ─────────────────────────────────────────────
# main() — mode x input source combinations
# ─────────────────────────────────────────────


class TestInputOptions:
    def run_manifest(self, tmp_path, extra_args):
        subprocess.run(
            [
                "python",
                "manifest_generator.py",
                "--outdir",
                str(tmp_path),
                "--name",
                "manifest.csv",
                *extra_args,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        return pd.read_csv(tmp_path / "manifest.csv")

    def assert_manifest(self, manifest, expected_columns, expected_ids):
        assert manifest.columns.tolist() == expected_columns
        assert len(manifest) == len(expected_ids)
        assert set(manifest["ID"]) == set(expected_ids)

    def test_short_from_dir(self, tmp_path):
        for name in ["sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz", "sampleB_R1.fastq.gz", "sampleB_R2.fastq.gz"]:
            (tmp_path / name).touch()

        manifest = self.run_manifest(tmp_path, ["--from_dir", str(tmp_path), "--mode", "short"])
        self.assert_manifest(manifest, ["ID", "R1", "R2"], ["sampleA", "sampleB"])

    def test_short_from_dir_recursive(self, tmp_path):
        sub = tmp_path / "subdir"
        sub.mkdir()
        for name in ["sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz", "sampleB_R1.fastq.gz", "sampleB_R2.fastq.gz"]:
            (sub / name).touch()
        for name in ["sampleC_R1.fastq.gz", "sampleC_R2.fastq.gz"]:
            (tmp_path / name).touch()

        manifest = self.run_manifest(
            tmp_path, ["--from_dir_recursive", str(tmp_path), "--max_depth", "1", "--mode", "short"]
        )
        self.assert_manifest(manifest, ["ID", "R1", "R2"], ["sampleA", "sampleB", "sampleC"])

    def test_long_from_dir(self, tmp_path):
        for name in ["sampleA.fastq.gz", "sampleB.fastq.gz"]:
            (tmp_path / name).touch()

        manifest = self.run_manifest(tmp_path, ["--from_dir", str(tmp_path), "--mode", "long"])
        self.assert_manifest(manifest, ["ID", "long_fastq", "genome_size"], ["sampleA", "sampleB"])

    def test_long_from_dir_recursive(self, tmp_path):
        sub = tmp_path / "subdir"
        sub.mkdir()
        for name in ["sampleA.fastq.gz", "sampleB.fastq.gz"]:
            (sub / name).touch()
        for name in ["sampleC.fastq.gz"]:
            (tmp_path / name).touch()

        manifest = self.run_manifest(
            tmp_path, ["--from_dir_recursive", str(tmp_path), "--max_depth", "1", "--mode", "long"]
        )
        self.assert_manifest(manifest, ["ID", "long_fastq", "genome_size"], ["sampleA", "sampleB", "sampleC"])

    def test_short_from_paths(self, tmp_path):
        for name in ["sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz", "sampleB_R1.fastq.gz", "sampleB_R2.fastq.gz"]:
            (tmp_path / name).touch()

        paths_file = tmp_path / "paths.txt"
        paths_file.write_text(
            "\n".join(
                str(tmp_path / name)
                for name in [
                    "sampleA_R1.fastq.gz",
                    "sampleA_R2.fastq.gz",
                    "sampleB_R1.fastq.gz",
                    "sampleB_R2.fastq.gz",
                ]
            )
        )

        manifest = self.run_manifest(tmp_path, ["--from_paths", str(paths_file), "--mode", "short"])
        self.assert_manifest(manifest, ["ID", "R1", "R2"], ["sampleA", "sampleB"])

    def test_long_from_paths(self, tmp_path):
        for name in ["sampleA.fastq.gz", "sampleB.fastq.gz"]:
            (tmp_path / name).touch()

        paths_file = tmp_path / "paths.txt"
        paths_file.write_text("sampleA.fastq.gz\n" "sampleB.fastq.gz")

        manifest = self.run_manifest(tmp_path, ["--from_paths", str(paths_file), "--mode", "long"])
        self.assert_manifest(manifest, ["ID", "long_fastq", "genome_size"], ["sampleA", "sampleB"])

    def test_hyrbid_from_paths_id(self, tmp_path):
        for name in [
            "sampleA_R1.fastq.gz",
            "sampleA_R2.fastq.gz",
            "sampleA_long.fastq.gz",
            "sampleB_R1.fastq.gz",
            "sampleB_R2.fastq.gz",
            "sampleB_long.fastq.gz",
        ]:
            (tmp_path / name).touch()

        paths_file = tmp_path / "paths.txt"
        paths_file.write_text(
            "ID,path\n"
            f"sampleA,{tmp_path}/sampleA_R1.fastq.gz\n"
            f"sampleA,{tmp_path}/sampleA_R2.fastq.gz\n"
            f"sampleA,{tmp_path}/sampleA_long.fastq.gz\n"
            f"sampleB,{tmp_path}/sampleB_R1.fastq.gz\n"
            f"sampleB,{tmp_path}/sampleB_R2.fastq.gz\n"
            f"sampleB,{tmp_path}/sampleB_long.fastq.gz\n"
        )

        manifest = self.run_manifest(tmp_path, ["--mode", "hybrid", "--from_paths_id", str(paths_file)])
        self.assert_manifest(manifest, ["ID", "R1", "R2", "long_fastq", "genome_size"], ["sampleA", "sampleB"])

    # TODO Need to mock the lgre tool to output a genome size (probably with a fake file fixture)
    # def test_genome_size_prediction(self, tmp_path):
    #     manifest = self.run_manifest(
    #         tmp_path, ["--mode", "long", "--from_dir", "testing/", "--genome_size", "--threads", "4"]
    #     )
    #     self.assert_manifest(manifest, ["ID", "long_fastq", "genome_size"], ["SRR38429962"])
    #     df = pd.read_csv(tmp_path / "manifest.csv")
    #     assert int(df.loc[0, "genome_size"]) >= 4400000

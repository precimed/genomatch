import json

import pytest

from test_prepare_variants import base_env, write_bim
from utils import read_tsv, run_py, run_py_with_env


@pytest.mark.parametrize(
    ("args", "message"),
    [
        (["--input", "source.bim", "--prefix", "work/@/prepared", "--output", "prepared"], "requires --input"),
        (["--input", "source.@.bim", "--prefix", "work/prepared", "--output", "prepared"], "requires --prefix"),
        (
            ["--input", "source.@.bim", "--prefix", "work/@/prepared", "--output", "prepared.@"],
            "requires non-sharded --output",
        ),
        (
            [
                "--input",
                "source.@.bim",
                "--prefix",
                "work/@/prepared",
                "--output",
                "prepared",
                "--dst-contig-naming",
                "plink_splitx",
            ],
            "does not support --dst-contig-naming=plink_splitx",
        ),
    ],
)
def test_prepare_variants_sharded_validates_sharded_contract(args, message):
    result = run_py("prepare_variants_sharded.py", *args, "--input-format", "bim")

    assert result.returncode != 0
    assert message in result.stderr


def test_prepare_variants_sharded_groups_x_tokens_and_concatenates_in_contig_rank_order(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "A" * 200, "chr2": "A" * 200, "chr10": "A" * 200, "chrX": "A" * 200},
        grch38_sequences={"chr1": "T" * 200, "chr2": "T" * 200, "chr10": "T" * 200, "chrX": "T" * 200},
    )
    source_template = tmp_path / "reference.@.bim"
    write_bim(tmp_path / "reference.10.bim", ["10\trs10\t0\t10\tG\tA"])
    write_bim(tmp_path / "reference.2.bim", ["2\trs2\t0\t2\tC\tA"])
    write_bim(tmp_path / "reference.NONPAR.bim", ["X\trs_nonpar\t0\t100\tG\tA"])
    write_bim(tmp_path / "reference.nonpar.bim", ["XY\trs_nonpar_lower\t0\t120\tC\tA"])
    prefix = tmp_path / "work" / "reference.@.prepared"
    output = tmp_path / "work" / "reference.prepared"

    result = run_py_with_env(
        "prepare_variants_sharded.py",
        env,
        "--input",
        source_template,
        "--input-format",
        "bim",
        "--prefix",
        prefix,
        "--output",
        output,
        "--shards",
        "10,2,NONPAR,nonpar",
        "--dst-build",
        "GRCh37",
        "--no-norm-indels",
    )

    assert result.returncode == 0, result.stderr
    assert "sort_variants.py" not in result.stderr

    command_lines = [line for line in result.stderr.splitlines() if "+ prepare_variants.py" in line]
    x_lines = [line for line in command_lines if str(tmp_path / "work" / "reference.X.prepared") in line]
    assert len(x_lines) == 1
    assert "--shards NONPAR,nonpar" in x_lines[0]
    assert "--contigs X" in x_lines[0]

    assert read_tsv(tmp_path / "work" / "reference.prepared.vmap") == [
        ["2", "2", "rs2", "C", "A", "2", "0", "identity"],
        ["10", "10", "rs10", "G", "A", "10", "0", "identity"],
        ["X", "100", "rs_nonpar", "G", "A", "NONPAR", "0", "identity"],
        ["X", "120", "rs_nonpar_lower", "C", "A", "nonpar", "0", "identity"],
    ]
    assert not (tmp_path / "work" / "reference.prepared.vmap.qc.tsv").exists()

    metadata = json.loads((tmp_path / "work" / "reference.prepared.vmap.meta.json").read_text(encoding="utf-8"))
    assert metadata["object_type"] == "variant_map"
    assert metadata["target"]["genome_build"] == "GRCh37"
    assert metadata["derived_from"] == str(source_template)


def test_prepare_variants_sharded_resume_noops_when_final_exists(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
        grch38_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
    )
    source_template = tmp_path / "reference.@.bim"
    write_bim(tmp_path / "reference.1.bim", ["1\trs1\t0\t1\tG\tA"])
    prefix = tmp_path / "work" / "reference.@.prepared"
    output = tmp_path / "work" / "reference.prepared"

    first = run_py_with_env(
        "prepare_variants_sharded.py",
        env,
        "--input",
        source_template,
        "--input-format",
        "bim",
        "--prefix",
        prefix,
        "--output",
        output,
        "--dst-build",
        "GRCh37",
    )
    assert first.returncode == 0, first.stderr

    second = run_py_with_env(
        "prepare_variants_sharded.py",
        env,
        "--input",
        source_template,
        "--input-format",
        "bim",
        "--prefix",
        prefix,
        "--output",
        output,
        "--dst-build",
        "GRCh37",
        "--resume",
    )

    assert second.returncode == 0, second.stderr
    assert second.stderr.strip().endswith(f"prepare_variants_sharded.py: wrote {tmp_path / 'work' / 'reference.prepared.vmap'}")
    assert "prepare_variants.py" not in second.stderr

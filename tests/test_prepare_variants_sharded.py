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
    ],
)
def test_prepare_variants_sharded_validates_sharded_contract(args, message):
    result = run_py("prepare_variants_sharded.py", *args, "--input-format", "bim")

    assert result.returncode != 0
    assert message in result.stderr


def test_prepare_variants_sharded_prepares_shards_and_sorts_final_target(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
        grch38_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
    )
    source_template = tmp_path / "reference.@.bim"
    write_bim(tmp_path / "reference.1.bim", ["1\trs1\t0\t1\tG\tA"])
    write_bim(tmp_path / "reference.2.bim", ["1\trs1_dup\t0\t1\tG\tA", "2\trs2\t0\t2\tC\tA"])
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
        "--dst-build",
        "GRCh37",
        "--no-norm-indels",
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "work" / "reference.1.prepared.vmap").exists()
    assert (tmp_path / "work" / "reference.2.prepared.vmap").exists()
    assert read_tsv(tmp_path / "work" / "reference.prepared.vmap") == [
        ["1", "1", "rs1", "G", "A", "1", "0", "identity"],
        ["2", "2", "rs2", "C", "A", "2", "1", "identity"],
    ]
    assert not (tmp_path / "work" / "reference.prepared.vmap.qc.tsv").exists()

    metadata = json.loads((tmp_path / "work" / "reference.prepared.vmap.meta.json").read_text(encoding="utf-8"))
    assert metadata["object_type"] == "variant_map"
    assert metadata["target"]["genome_build"] == "GRCh37"
    assert metadata["derived_from"] == str(source_template)

    sort_prefix = tmp_path / "work" / "reference.all_targets.prepared.sort_tmp"
    assert "sort_variants.py" in result.stderr
    assert "--drop-duplicates" in result.stderr
    assert str(sort_prefix) in result.stderr
    assert "dropped 1 duplicate target rows" in result.stderr


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
    assert "sort_variants.py" not in second.stderr

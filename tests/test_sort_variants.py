import json
import os

import pytest

from utils import read_tsv, run_py, run_py_with_env, write_json, write_lines


def test_sort_variants_vtable_sorts_by_declared_coordinate_order(tmp_path):
    source = tmp_path / "source.vtable"
    out = tmp_path / "sorted.vtable"
    write_lines(
        source,
        [
            "X\t7\trx\tA\tC",
            "1\t200\tr2\tC\tT",
            "1\t100\tr1\tA\tG",
            "MT\t3\trm\tA\tT",
            "1\t100\tr1b\tA\tC",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("sort_variants.py", "--input", source, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "100", "r1", "A", "G"],
        ["1", "100", "r1b", "A", "C"],
        ["1", "200", "r2", "C", "T"],
        ["X", "7", "rx", "A", "C"],
        ["MT", "3", "rm", "A", "T"],
    ]
    assert json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8")) == {
        "object_type": "variant_table",
        "genome_build": "GRCh37",
        "contig_naming": "ncbi",
    }


def test_sort_variants_vmap_sorts_target_rows_and_preserves_provenance(tmp_path):
    source = tmp_path / "source.vmap"
    out = tmp_path / "sorted.vmap"
    write_lines(
        source,
        [
            "2\t5\tr2\tA\tG\tsh2\t9\tswap",
            "1\t4\tr1b\tC\tT\tsh1\t8\tflip",
            "1\t4\tr1a\tA\tG\tsh1\t7\tidentity",
            "X\t1\trx\tA\tC\tsh3\t3\tflip_swap",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh38", "contig_naming": "ncbi"}},
    )

    result = run_py("sort_variants.py", "--input", source, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "4", "r1b", "C", "T", "sh1", "8", "flip"],
        ["1", "4", "r1a", "A", "G", "sh1", "7", "identity"],
        ["2", "5", "r2", "A", "G", "sh2", "9", "swap"],
        ["X", "1", "rx", "A", "C", "sh3", "3", "flip_swap"],
    ]
    assert json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8")) == {
        "object_type": "variant_map",
        "target": {"genome_build": "GRCh38", "contig_naming": "ncbi"},
    }


def test_sort_variants_removes_stale_qc_sidecar(tmp_path):
    source = tmp_path / "source.vtable"
    out = tmp_path / "sorted.vtable"
    write_lines(source, ["1\t2\tr2\tC\tT", "1\t1\tr1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )
    write_lines(out.with_name(out.name + ".qc.tsv"), ["stale"])

    result = run_py("sort_variants.py", "--input", source, "--output", out)

    assert result.returncode == 0, result.stderr
    assert not out.with_name(out.name + ".qc.tsv").exists()


def test_sort_variants_drop_duplicates_handles_non_adjacent_identity_with_same_coordinate(tmp_path):
    source = tmp_path / "source.vmap"
    out = tmp_path / "sorted.vmap"
    write_lines(
        source,
        [
            "1\t10\tfirst\tA\tG\tsh1\t1\tidentity",
            "1\t10\tother\tC\tT\tsh2\t2\tidentity",
            "1\t10\tduplicate\tA\tG\tsh3\t3\tidentity",
            "1\t11\tlater\tA\tC\tsh4\t4\tidentity",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh38", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "sort_variants.py",
        {"GENOMATCH_SORT_CHUNK_LINES": "2"},
        "--input",
        source,
        "--output",
        out,
        "--drop-duplicates",
        "--prefix",
        tmp_path / "scratch" / "sort_tmp",
    )

    assert result.returncode == 0, result.stderr
    assert "dropped 1 duplicate target rows" in result.stderr
    assert read_tsv(out) == [
        ["1", "10", "first", "A", "G", "sh1", "1", "identity"],
        ["1", "10", "other", "C", "T", "sh2", "2", "identity"],
        ["1", "11", "later", "A", "C", "sh4", "4", "identity"],
    ]


def test_sort_variants_rejects_invalid_chunk_env(tmp_path):
    source = tmp_path / "source.vtable"
    out = tmp_path / "sorted.vtable"
    write_lines(source, ["1\t1\tr1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "sort_variants.py",
        {"GENOMATCH_SORT_CHUNK_LINES": "0"},
        "--input",
        source,
        "--output",
        out,
    )

    assert result.returncode != 0
    assert "GENOMATCH_SORT_CHUNK_LINES must be a positive integer" in result.stderr


def test_sort_variants_unwritable_scratch_prefix_fails_clearly(tmp_path):
    source = tmp_path / "source.vtable"
    out = tmp_path / "sorted.vtable"
    scratch_parent = tmp_path / "scratch"
    scratch_parent.mkdir()
    write_lines(source, ["1\t2\tr2\tC\tT", "1\t1\tr1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    scratch_parent.chmod(0o500)
    try:
        if os.access(scratch_parent, os.W_OK):
            pytest.skip("scratch directory remains writable in this environment")
        result = run_py(
            "sort_variants.py",
            "--input",
            source,
            "--output",
            out,
            "--prefix",
            scratch_parent / "sort_tmp",
        )
    finally:
        scratch_parent.chmod(0o700)

    assert result.returncode != 0
    assert "Permission denied" in result.stderr
    assert not out.exists()


def test_sort_variants_external_sort_merges_many_runs(tmp_path):
    source = tmp_path / "source.vtable"
    out = tmp_path / "sorted.vtable"
    write_lines(source, [f"1\t{idx}\tr{idx}\tA\tG" for idx in range(130, 0, -1)])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "sort_variants.py",
        {"GENOMATCH_SORT_CHUNK_LINES": "1"},
        "--input",
        source,
        "--output",
        out,
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [[str(1), str(idx), f"r{idx}", "A", "G"] for idx in range(1, 131)]


def test_sort_variants_default_scratch_is_output_adjacent(tmp_path):
    source = tmp_path / "source.vtable"
    out = tmp_path / "nested" / "sorted.vtable"
    write_lines(source, ["1\t2\tr2\tC\tT", "1\t1\tr1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "sort_variants.py",
        {"GENOMATCH_SORT_CHUNK_LINES": "1"},
        "--input",
        source,
        "--output",
        out,
    )

    assert result.returncode == 0, result.stderr
    assert str(out.with_name(out.name + ".sort_tmp.")) in result.stderr

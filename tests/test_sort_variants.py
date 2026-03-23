import json

from utils import read_tsv, run_py, write_json, write_lines


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

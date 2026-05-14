import json

from utils import read_tsv, run_py, write_json, write_lines


def test_union_variants_deduplicates_target_keys_assigns_target_ids_and_sorts(tmp_path):
    first = tmp_path / "a.vmap"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(
        first,
        [
            "X\t7\trx\tA\tC\tshx\t0\tidentity",
            "1\t200\tlate\tC\tT\tsh1\t1\tidentity",
            "1\t100\tfirst-pos\tA\tG\tsh1\t2\tidentity",
        ],
    )
    write_lines(
        second,
        [
            "1\t100\tsecond-copy\tA\tG",
            "1\t100\tsecond-copy-dup\tA\tG",
            "1\t100\tother-alleles\tA\tC",
            "MT\t3\tmt\tA\tT",
        ],
    )
    write_json(
        first.with_name(first.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )
    write_json(
        second.with_name(second.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("union_variants.py", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "100", "1:100:A:C", "A", "C"],
        ["1", "100", "1:100:A:G", "A", "G"],
        ["1", "200", "1:200:C:T", "C", "T"],
        ["X", "7", "X:7:A:C", "A", "C"],
        ["MT", "3", "MT:3:A:T", "A", "T"],
    ]
    assert json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))["variants_count"] == 5
    assert f"loaded 3 variants from {first}" in result.stderr
    assert f"after unioning {first}, 3 variants accumulated" in result.stderr
    assert f"loaded 4 variants from {second}" in result.stderr
    assert f"after unioning {second}, 7 variants accumulated" in result.stderr
    assert "variants_count missing from metadata" in result.stderr
    assert "using CLI input order" in result.stderr


def test_union_variants_rewrites_ids_for_non_overlapping_inputs(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(first, ["1\t100\tinput-a\tA\tG"])
    write_lines(second, ["1\t200\tinput-b\tC\tT"])
    write_json(first.with_name(first.name + ".meta.json"), meta)
    write_json(second.with_name(second.name + ".meta.json"), meta)

    result = run_py("union_variants.py", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "100", "1:100:A:G", "A", "G"],
        ["1", "200", "1:200:C:T", "C", "T"],
    ]


def test_union_variants_uses_count_order_when_all_counts_are_present(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t300\tinput-a\tA\tC", "1\t100\tinput-b\tA\tG"])
    write_lines(second, ["1\t200\tinput-c\tC\tT"])
    write_json(
        first.with_name(first.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi", "variants_count": 2},
    )
    write_json(
        second.with_name(second.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi", "variants_count": 1},
    )

    result = run_py("union_variants.py", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "100", "1:100:A:G", "A", "G"],
        ["1", "200", "1:200:C:T", "C", "T"],
        ["1", "300", "1:300:A:C", "A", "C"],
    ]
    assert result.stderr.index(f"loaded 1 variants from {second}") < result.stderr.index(f"loaded 2 variants from {first}")
    assert "variants_count missing from metadata" not in result.stderr


def test_union_variants_requires_at_least_two_inputs(tmp_path):
    source = tmp_path / "a.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["1\t100\tid1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("union_variants.py", source, "--output", out)

    assert result.returncode != 0
    assert "at least two inputs" in result.stderr


def test_union_variants_rejects_legacy_inputs_flag(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["1\t200\tid2\tC\tT"])
    write_json(first.with_name(first.name + ".meta.json"), meta)
    write_json(second.with_name(second.name + ".meta.json"), meta)

    result = run_py("union_variants.py", "--inputs", first, second, "--output", out)

    assert result.returncode != 0
    assert "unrecognized arguments" in result.stderr


def test_union_variants_checks_metadata_before_loading_rows(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["not\ta\tvalid\tvtable\trow"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": "ncbi"})

    result = run_py("union_variants.py", first, second, "--output", out)

    assert result.returncode != 0
    assert "same genome_build" in result.stderr
    assert f"{first}: genome_build=GRCh37, contig_naming=ncbi" in result.stderr
    assert f"{second}: genome_build=GRCh38, contig_naming=ncbi" in result.stderr
    assert "invalid vtable row" not in result.stderr


def test_union_variants_checks_contig_naming_before_loading_rows(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["not\ta\tvalid\tvtable\trow"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ucsc"})

    result = run_py("union_variants.py", first, second, "--output", out)

    assert result.returncode != 0
    assert "same contig_naming" in result.stderr
    assert f"{first}: genome_build=GRCh37, contig_naming=ncbi" in result.stderr
    assert f"{second}: genome_build=GRCh37, contig_naming=ucsc" in result.stderr
    assert "invalid vtable row" not in result.stderr


def test_union_variants_fails_cleanly_on_missing_contig_naming(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["not\ta\tvalid\tvtable\trow"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py("union_variants.py", first, second, "--output", out)

    assert result.returncode != 0
    assert "contig_naming" in result.stderr
    assert "normalize_contigs.py" in result.stderr
    assert f"{first}: genome_build=GRCh37, contig_naming=<missing>" in result.stderr
    assert "invalid vtable row" not in result.stderr


def test_union_variants_rejects_template_paths(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["1\t200\tid2\tC\tT"])
    write_json(first.with_name(first.name + ".meta.json"), meta)
    write_json(second.with_name(second.name + ".meta.json"), meta)

    result = run_py("union_variants.py", first, "b.@.vtable", "--output", out)

    assert result.returncode != 0
    assert "does not accept '@' paths" in result.stderr

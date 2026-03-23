from utils import read_tsv, run_py, write_json, write_lines


def test_union_variants_deduplicates_first_occurrence_and_sorts(tmp_path):
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

    result = run_py("union_variants.py", "--inputs", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "100", "first-pos", "A", "G"],
        ["1", "100", "other-alleles", "A", "C"],
        ["1", "200", "late", "C", "T"],
        ["X", "7", "rx", "A", "C"],
        ["MT", "3", "mt", "A", "T"],
    ]
    assert f"loaded 3 variants from {first}" in result.stderr
    assert f"after unioning {first}, 3 variants accumulated" in result.stderr
    assert f"loaded 4 variants from {second}" in result.stderr
    assert f"after unioning {second}, 5 variants accumulated" in result.stderr


def test_union_variants_requires_at_least_two_inputs(tmp_path):
    source = tmp_path / "a.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["1\t100\tid1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("union_variants.py", "--inputs", source, "--output", out)

    assert result.returncode != 0
    assert "at least two inputs" in result.stderr


def test_union_variants_checks_metadata_before_loading_rows(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["not\ta\tvalid\tvtable\trow"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": "ncbi"})

    result = run_py("union_variants.py", "--inputs", first, second, "--output", out)

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

    result = run_py("union_variants.py", "--inputs", first, second, "--output", out)

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

    result = run_py("union_variants.py", "--inputs", first, second, "--output", out)

    assert result.returncode != 0
    assert "contig_naming" in result.stderr
    assert "normalize_contigs.py" in result.stderr
    assert f"{first}: genome_build=GRCh37, contig_naming=<missing>" in result.stderr
    assert "invalid vtable row" not in result.stderr

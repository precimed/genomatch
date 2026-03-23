from utils import read_tsv, run_py, write_json, write_lines


def test_intersect_variants_preserves_first_input_order(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(first, ["1\t200\tid2\tC\tT", "1\t100\tid1\tA\tG", "1\t300\tid3\tA\tC"])
    write_lines(second, ["1\t100\tx\tA\tG", "1\t200\ty\tC\tT"])
    write_json(first.with_name(first.name + ".meta.json"), meta)
    write_json(second.with_name(second.name + ".meta.json"), meta)
    result = run_py("intersect_variants.py", "--inputs", first, second, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "200", "id2", "C", "T"], ["1", "100", "id1", "A", "G"]]
    assert f"loaded 3 variants from {first}" in result.stderr
    assert f"loaded 2 variants from {second}" in result.stderr
    assert f"after intersecting {second}, 2 variants remain" in result.stderr


def test_intersect_variants_checks_metadata_before_loading_rows(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["not\ta\tvalid\tvtable\trow"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": "ncbi"})

    result = run_py("intersect_variants.py", "--inputs", first, second, "--output", out)

    assert result.returncode != 0
    assert "same genome_build" in result.stderr
    assert f"{first}: genome_build=GRCh37, contig_naming=ncbi" in result.stderr
    assert f"{second}: genome_build=GRCh38, contig_naming=ncbi" in result.stderr
    assert "invalid vtable row" not in result.stderr


def test_intersect_variants_checks_contig_naming_before_loading_rows(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["not\ta\tvalid\tvtable\trow"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ucsc"})

    result = run_py("intersect_variants.py", "--inputs", first, second, "--output", out)

    assert result.returncode != 0
    assert "same contig_naming" in result.stderr
    assert f"{first}: genome_build=GRCh37, contig_naming=ncbi" in result.stderr
    assert f"{second}: genome_build=GRCh37, contig_naming=ucsc" in result.stderr
    assert "invalid vtable row" not in result.stderr


def test_intersect_variants_fails_cleanly_on_missing_contig_naming(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["not\ta\tvalid\tvtable\trow"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py("intersect_variants.py", "--inputs", first, second, "--output", out)

    assert result.returncode != 0
    assert "contig_naming" in result.stderr
    assert "normalize_contigs.py" in result.stderr
    assert f"{first}: genome_build=GRCh37, contig_naming=<missing>" in result.stderr
    assert "invalid vtable row" not in result.stderr

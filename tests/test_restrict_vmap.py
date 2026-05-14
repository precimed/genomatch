from utils import read_tsv, run_py, write_json, write_lines


def vmap_meta() -> dict:
    return {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}}


def vtable_meta() -> dict:
    return {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}


def write_vmap(path, rows, metadata=None):
    write_lines(path, rows)
    write_json(path.with_name(path.name + ".meta.json"), metadata or vmap_meta())


def write_vtable(path, rows, metadata=None):
    write_lines(path, rows)
    write_json(path.with_name(path.name + ".meta.json"), metadata or vtable_meta())


def test_restrict_vmap_preserves_source_order_ids_provenance_and_allele_op(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(
        source,
        [
            "1\t300\tsrc3\tA\tC\t.\t2\tflip",
            "1\t100\tsrc1\tA\tG\t.\t0\tidentity",
            "1\t200\tsrc2\tC\tT\t.\t1\tswap",
        ],
    )
    write_vtable(restriction, ["1\t100\tother1\tA\tG", "1\t300\tother3\tA\tC"])

    result = run_py("restrict_vmap.py", source, restriction, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "300", "src3", "A", "C", ".", "2", "flip"],
        ["1", "100", "src1", "A", "G", ".", "0", "identity"],
    ]
    assert read_tsv(source) == [
        ["1", "300", "src3", "A", "C", ".", "2", "flip"],
        ["1", "100", "src1", "A", "G", ".", "0", "identity"],
        ["1", "200", "src2", "C", "T", ".", "1", "swap"],
    ]
    assert "variants_count missing from metadata" in result.stderr
    assert "using CLI input order" in result.stderr


def test_restrict_vmap_uses_count_order_but_preserves_source_output_order(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(
        source,
        [
            "1\t300\tsrc3\tA\tC\t.\t2\tflip",
            "1\t100\tsrc1\tA\tG\t.\t0\tidentity",
            "1\t200\tsrc2\tC\tT\t.\t1\tswap",
        ],
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}, "variants_count": 3},
    )
    write_vtable(
        restriction,
        ["1\t100\tother1\tA\tG"],
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi", "variants_count": 1},
    )

    result = run_py("restrict_vmap.py", source, restriction, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "src1", "A", "G", ".", "0", "identity"]]
    assert result.stderr.index(f"loaded 1 restriction variants from {restriction}") < result.stderr.index(
        f"loaded 3 source variants from {source}"
    )
    assert "variants_count missing from metadata" not in result.stderr


def test_restrict_vmap_ignores_restriction_vmap_provenance(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vmap"
    out = tmp_path / "out.vmap"
    write_vmap(
        source,
        [
            "1\t100\tsource_id\tA\tG\tsource_shard\t7\tidentity",
            "1\t200\tsource_only\tC\tT\tsource_shard\t8\tidentity",
        ],
    )
    write_vmap(restriction, ["1\t100\trestriction_id\tA\tG\tignored_shard\t99\tflip_swap"])

    result = run_py("restrict_vmap.py", source, restriction, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "source_id", "A", "G", "source_shard", "7", "identity"]]


def test_restrict_vmap_multiple_restrictions_are_intersection_not_union(tmp_path):
    source = tmp_path / "source.vmap"
    first = tmp_path / "first.vtable"
    second = tmp_path / "second.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(
        source,
        [
            "1\t100\tboth\tA\tG\t.\t0\tidentity",
            "1\t200\tfirst_only\tC\tT\t.\t1\tidentity",
            "1\t300\tsecond_only\tA\tC\t.\t2\tidentity",
        ],
    )
    write_vtable(first, ["1\t100\tx\tA\tG", "1\t200\ty\tC\tT"])
    write_vtable(second, ["1\t100\tx\tA\tG", "1\t300\tz\tA\tC"])

    result = run_py("restrict_vmap.py", source, first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "both", "A", "G", ".", "0", "identity"]]


def test_restrict_vmap_source_duplicate_target_keys_fail_before_restriction(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(
        source,
        [
            "1\t100\tfirst\tA\tG\t.\t0\tidentity",
            "1\t100\tsecond\tA\tG\t.\t1\tidentity",
        ],
    )
    write_vtable(restriction, ["1\t100\tx\tA\tG"])

    result = run_py("restrict_vmap.py", source, restriction, "--output", out)

    assert result.returncode != 0
    assert "duplicate chrom:pos:a1:a2 in vmap target rows" in result.stderr


def test_restrict_vmap_restriction_vmap_duplicate_target_keys_fail(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vmap"
    out = tmp_path / "out.vmap"
    write_vmap(source, ["1\t100\tid1\tA\tG\t.\t0\tidentity"])
    write_vmap(
        restriction,
        [
            "1\t100\tfirst\tA\tG\t.\t0\tidentity",
            "1\t100\tsecond\tA\tG\t.\t1\tidentity",
        ],
    )

    result = run_py("restrict_vmap.py", source, restriction, "--output", out)

    assert result.returncode != 0
    assert "duplicate chrom:pos:a1:a2 in vmap target rows" in result.stderr


def test_restrict_vmap_allows_duplicate_vtable_membership_keys(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(source, ["1\t100\tid1\tA\tG\t.\t0\tidentity"])
    write_vtable(restriction, ["1\t100\tfirst\tA\tG", "1\t100\tsecond\tA\tG"])

    result = run_py("restrict_vmap.py", source, restriction, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "id1", "A", "G", ".", "0", "identity"]]


def test_restrict_vmap_metadata_mismatch_fails_clearly(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(source, ["1\t100\tid1\tA\tG\t.\t0\tidentity"])
    write_vtable(
        restriction,
        ["1\t100\tid1\tA\tG"],
        {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": "ncbi"},
    )

    result = run_py("restrict_vmap.py", source, restriction, "--output", out)

    assert result.returncode != 0
    assert "same genome_build" in result.stderr
    assert f"{source}: genome_build=GRCh37, contig_naming=ncbi" in result.stderr
    assert f"{restriction}: genome_build=GRCh38, contig_naming=ncbi" in result.stderr


def test_restrict_vmap_missing_contig_naming_fails_clearly(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(
        source,
        ["1\t100\tid1\tA\tG\t.\t0\tidentity"],
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37"}},
    )
    write_vtable(restriction, ["1\t100\tid1\tA\tG"])

    result = run_py("restrict_vmap.py", source, restriction, "--output", out)

    assert result.returncode != 0
    assert "contig_naming" in result.stderr
    assert "normalize_contigs.py" in result.stderr
    assert f"{source}: genome_build=GRCh37, contig_naming=<missing>" in result.stderr


def test_restrict_vmap_supports_non_snp_biallelic_exact_restriction(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(
        source,
        [
            "1\t100\tindel\tAC\tA\t.\t0\tidentity",
            "1\t200\tsnp\tG\tT\t.\t1\tidentity",
        ],
    )
    write_vtable(restriction, ["1\t100\tother\tAC\tA"])

    result = run_py("restrict_vmap.py", source, restriction, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "indel", "AC", "A", ".", "0", "identity"]]


def test_restrict_vmap_rejects_missing_provenance_source_rows(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(source, ["1\t100\tmissing\tA\tG\t.\t-1\tmissing"])
    write_vtable(restriction, ["1\t100\tid1\tA\tG"])

    result = run_py("restrict_vmap.py", source, restriction, "--output", out)

    assert result.returncode != 0
    assert "mapped-only .vmap" in result.stderr


def test_restrict_vmap_rejects_template_paths_and_sort_flag(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(source, ["1\t100\tid1\tA\tG\t.\t0\tidentity"])
    write_vtable(restriction, ["1\t100\tid1\tA\tG"])

    template_result = run_py("restrict_vmap.py", source, "membership.@.vtable", "--output", out)
    sort_result = run_py("restrict_vmap.py", source, restriction, "--output", out, "--sort")

    assert template_result.returncode != 0
    assert "does not accept '@' paths" in template_result.stderr
    assert sort_result.returncode != 0
    assert "unrecognized arguments: --sort" in sort_result.stderr


def test_restrict_vmap_rejects_legacy_source_and_targets_flags(tmp_path):
    source = tmp_path / "source.vmap"
    restriction = tmp_path / "membership.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(source, ["1\t100\tid1\tA\tG\t.\t0\tidentity"])
    write_vtable(restriction, ["1\t100\tid1\tA\tG"])

    result = run_py("restrict_vmap.py", "--source", source, "--targets", restriction, "--output", out)

    assert result.returncode != 0
    assert "unrecognized arguments" in result.stderr

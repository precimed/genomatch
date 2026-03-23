from utils import read_tsv, run_py, write_json, write_lines


def test_match_vmap_to_target_identity_swap_and_missing(tmp_path):
    source = tmp_path / "src.vmap"
    target = tmp_path / "tgt.vtable"
    out = tmp_path / "map.vmap"
    write_lines(source, ["1\t100\trs1\tA\tG\t.\t0\tidentity", "1\t200\trs2\tC\tT\t.\t1\tidentity"])
    write_lines(target, ["1\t200\tt1\tT\tC", "1\t100\tt2\tA\tG", "1\t300\tt3\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )
    write_json(target.with_name(target.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    result = run_py("match_vmap_to_target.py", "--source", source, "--target", target, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "200", "t1", "T", "C", ".", "1", "swap"],
        ["1", "100", "t2", "A", "G", ".", "0", "identity"],
        ["1", "300", "t3", "A", "C", ".", "-1", "missing"],
    ]


def test_match_vmap_to_target_supports_non_snv_identity_and_swap(tmp_path):
    source = tmp_path / "src.vmap"
    target = tmp_path / "tgt.vtable"
    out = tmp_path / "map.vmap"
    meta_vmap = {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}}
    meta_vtable = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(source, ["1\t100\trs1\tA\tAT\t.\t0\tidentity", "1\t200\trs2\tCG\tC\t.\t1\tidentity"])
    write_lines(target, ["1\t100\tt1\tA\tAT", "1\t200\tt2\tC\tCG"])
    write_json(source.with_name(source.name + ".meta.json"), meta_vmap)
    write_json(target.with_name(target.name + ".meta.json"), meta_vtable)
    result = run_py("match_vmap_to_target.py", "--source", source, "--target", target, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "100", "t1", "A", "AT", ".", "0", "identity"],
        ["1", "200", "t2", "C", "CG", ".", "1", "swap"],
    ]


def test_match_vmap_to_target_requires_source_vmap(tmp_path):
    source = tmp_path / "src.vtable"
    target = tmp_path / "tgt.vmap"
    out = tmp_path / "map.vmap"
    write_lines(source, ["1\t100\trs1\tA\tG"])
    write_lines(target, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(source.with_name(source.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(target.with_name(target.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("match_vmap_to_target.py", "--source", source, "--target", target, "--output", out)
    assert result.returncode != 0
    assert "--source" in result.stderr or ".vmap" in result.stderr


def test_match_vmap_to_target_accepts_target_vmap_and_warns_that_provenance_is_ignored(tmp_path):
    source = tmp_path / "src.vmap"
    target = tmp_path / "tgt.vmap"
    out = tmp_path / "map.vmap"
    write_lines(source, ["1\t100\trs_source\tA\tG\tshard1\t8\tswap", "1\t200\trs2\tC\tT\tshard2\t9\tidentity"])
    write_lines(
        target,
        [
            "1\t200\tt2\tT\tC\ttargetShard\t55\tflip",
            "1\t100\tt1\tG\tA\totherShard\t77\tidentity",
            "1\t300\tt3\tA\tC\t.\t-1\tmissing",
        ],
    )
    write_json(source.with_name(source.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})
    write_json(target.with_name(target.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("match_vmap_to_target.py", "--source", source, "--target", target, "--output", out)

    assert result.returncode == 0, result.stderr
    assert "provenance is ignored" in result.stderr
    assert "convert_vmap_to_target.py" in result.stderr
    assert read_tsv(out) == [
        ["1", "200", "t2", "T", "C", "shard2", "9", "swap"],
        ["1", "100", "t1", "G", "A", "shard1", "8", "identity"],
        ["1", "300", "t3", "A", "C", ".", "-1", "missing"],
    ]


def test_match_vmap_to_target_rejects_mismatched_contig_naming(tmp_path):
    source = tmp_path / "src.vmap"
    target = tmp_path / "tgt.vtable"
    out = tmp_path / "map.vmap"
    write_lines(source, ["1\t100\trs1\tA\tG\t.\t0\tidentity"])
    write_lines(target, ["chr1\t100\tt1\tA\tG"])
    write_json(source.with_name(source.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})
    write_json(target.with_name(target.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ucsc"})

    result = run_py("match_vmap_to_target.py", "--source", source, "--target", target, "--output", out)
    assert result.returncode != 0
    assert "contig_naming" in result.stderr


def test_match_vmap_to_target_rejects_duplicate_target_rows(tmp_path):
    source = tmp_path / "src.vmap"
    target = tmp_path / "tgt.vtable"
    out = tmp_path / "map.vmap"
    meta_vmap = {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}}
    meta_vtable = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(source, ["1\t100\trs1\tA\tG\t.\t0\tidentity"])
    write_lines(target, ["1\t100\tt1\tA\tG", "1\t100\tt2\tA\tG"])
    write_json(source.with_name(source.name + ".meta.json"), meta_vmap)
    write_json(target.with_name(target.name + ".meta.json"), meta_vtable)

    result = run_py("match_vmap_to_target.py", "--source", source, "--target", target, "--output", out)
    assert result.returncode != 0
    assert "duplicate chrom:pos:a1:a2" in result.stderr.lower()


def test_match_vmap_to_target_rejects_unknown_target_contig(tmp_path):
    source = tmp_path / "src.vmap"
    target = tmp_path / "tgt.vtable"
    out = tmp_path / "map.vmap"
    meta_vmap = {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}}
    meta_vtable = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(source, ["1\t100\trs1\tA\tG\t.\t0\tidentity"])
    write_lines(target, ["unknown\t100\tt1\tA\tG"])
    write_json(source.with_name(source.name + ".meta.json"), meta_vmap)
    write_json(target.with_name(target.name + ".meta.json"), meta_vtable)

    result = run_py("match_vmap_to_target.py", "--source", source, "--target", target, "--output", out)
    assert result.returncode != 0
    assert "normalize_contigs.py" in result.stderr


def test_match_vmap_to_target_rejects_rows_inconsistent_with_declared_contig_naming(tmp_path):
    source = tmp_path / "src.vmap"
    target = tmp_path / "tgt.vtable"
    out = tmp_path / "map.vmap"
    write_lines(source, ["chr1\t100\trs1\tA\tG\t.\t0\tidentity"])
    write_lines(target, ["1\t100\tt1\tA\tG"])
    write_json(source.with_name(source.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})
    write_json(target.with_name(target.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py("match_vmap_to_target.py", "--source", source, "--target", target, "--output", out)
    assert result.returncode != 0
    assert "normalize_contigs.py" in result.stderr


def test_match_vmap_to_target_ignores_rsid_for_matching_and_preserves_provenance(tmp_path):
    source = tmp_path / "src.vmap"
    target = tmp_path / "tgt.vtable"
    out = tmp_path / "map.vmap"
    write_lines(source, ["1\t100\trs_source\tA\tG\tshard1\t8\tswap"])
    write_lines(target, ["1\t100\trs_target\tG\tA"])
    write_json(source.with_name(source.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})
    write_json(target.with_name(target.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py("match_vmap_to_target.py", "--source", source, "--target", target, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "rs_target", "G", "A", "shard1", "8", "identity"]]


def test_match_vmap_to_target_uses_first_source_duplicate_match(tmp_path):
    source = tmp_path / "src.vmap"
    target = tmp_path / "tgt.vtable"
    out = tmp_path / "map.vmap"
    write_lines(
        source,
        [
            "1\t100\trs_first\tA\tG\tshard1\t8\tidentity",
            "1\t100\trs_second\tG\tA\tshard2\t9\tidentity",
        ],
    )
    write_lines(target, ["1\t100\tt1\tA\tG"])
    write_json(source.with_name(source.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})
    write_json(target.with_name(target.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py("match_vmap_to_target.py", "--source", source, "--target", target, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "t1", "A", "G", "shard1", "8", "identity"]]


def test_match_vmap_to_target_rejects_template_paths(tmp_path):
    source_template = tmp_path / "src.@.vmap"
    target = tmp_path / "tgt.vtable"
    output_template = tmp_path / "out.@.vmap"
    source1 = tmp_path / "src.1.vmap"
    write_lines(source1, ["1\t100\trs1\tA\tG\t.\t0\tidentity"])
    write_json(source1.with_name(source1.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})
    write_lines(target, ["1\t100\tt1\tA\tG"])
    write_json(target.with_name(target.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py(
        "match_vmap_to_target.py",
        "--source",
        source_template,
        "--target",
        target,
        "--output",
        output_template,
    )
    assert result.returncode != 0
    assert "does not accept '@' paths" in result.stderr

import os

from genomatch.apply_vmap_bfile import has_identity_sample_axis_remap
from utils import read_bed_matrix, read_bim, read_fam_sample_count, run_py, run_py_with_env, write_bfile, write_json, write_lines


def test_has_identity_sample_axis_remap_requires_full_identity_mapping():
    assert has_identity_sample_axis_remap((0, 1, 2), 3) is True
    assert has_identity_sample_axis_remap((1, 2, 0), 3) is False
    assert has_identity_sample_axis_remap((0, 1), 3) is False
    assert has_identity_sample_axis_remap((0, -1, 2), 3) is False


def test_apply_vmap_to_bfile_preserves_target_order_and_missing_fill_for_single_file_provenance(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(
        source_prefix,
        ["1\trs1\t0\t100\tA\tG", "1\trs2\t0\t200\tC\tT", "1\trs3\t0\t300\tG\tA"],
        ["S1 S1 0 0 0 -9", "S2 S2 0 0 0 -9"],
        [[0, 3], [0, 2], [3, 0]],
    )
    write_lines(
        vmap,
        [
            "1\t200\tt2\tT\tC\t.\t1\tswap",
            "1\t100\tt1\tA\tG\t.\t0\tidentity",
            "1\t400\tt3\tA\tC\t.\t-1\tmissing",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})
    os.utime(source_prefix.with_suffix(".fam"), (1, 1))

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert [(row.chrom, row.snp, row.cm, row.bp, row.a1, row.a2) for row in read_bim(out_prefix.with_suffix(".bim"))] == [
        ("1", "1:200:T:C", "0", "200", "T", "C"),
        ("1", "1:100:A:G", "0", "100", "A", "G"),
        ("1", "1:400:A:C", "0", "400", "A", "C"),
    ]
    assert out_prefix.with_suffix(".fam").read_text(encoding="utf-8") == source_prefix.with_suffix(".fam").read_text(encoding="utf-8")
    assert out_prefix.with_suffix(".fam").stat().st_mtime > source_prefix.with_suffix(".fam").stat().st_mtime
    matrix = read_bed_matrix(out_prefix.with_suffix(".bed"), 2, 3)
    assert matrix == [[3, 2], [0, 3], [1, 1]]


def test_apply_vmap_to_bfile_only_mapped_target_drops_unmatched_rows(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(
        source_prefix,
        ["1\trs1\t0\t100\tA\tG", "1\trs2\t0\t200\tC\tT", "1\trs3\t0\t300\tG\tA"],
        ["S1 S1 0 0 0 -9", "S2 S2 0 0 0 -9"],
        [[0, 3], [0, 2], [3, 0]],
    )
    write_lines(
        vmap,
        [
            "1\t200\tt2\tT\tC\t.\t1\tswap",
            "1\t400\tt_missing\tA\tC\t.\t-1\tmissing",
            "1\t100\tt1\tA\tG\t.\t0\tidentity",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py(
        "apply_vmap_to_bfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--only-mapped-target",
    )
    assert result.returncode == 0, result.stderr
    assert [(row.chrom, row.snp, row.cm, row.bp, row.a1, row.a2) for row in read_bim(out_prefix.with_suffix(".bim"))] == [
        ("1", "1:200:T:C", "0", "200", "T", "C"),
        ("1", "1:100:A:G", "0", "100", "A", "G"),
    ]
    matrix = read_bed_matrix(out_prefix.with_suffix(".bed"), 2, 2)
    assert matrix == [[3, 2], [0, 3]]


def test_apply_vmap_to_bfile_retain_snp_id_uses_vmap_id(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(
        source_prefix,
        ["1\trs_source\t0\t100\tA\tG"],
        ["S1 S1 0 0 0 -9"],
        [[0]],
    )
    write_lines(vmap, ["1\t100\ttarget_id\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py(
        "apply_vmap_to_bfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--retain-snp-id",
    )

    assert result.returncode == 0, result.stderr
    assert [row.snp for row in read_bim(out_prefix.with_suffix(".bim"))] == ["target_id"]


def test_apply_vmap_to_bfile_treats_flip_swap_as_genotype_swap(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 0 -9", "S2 S2 0 0 0 -9"], [[0, 3]])
    write_lines(vmap, ["1\t100\tt1\tC\tT\t.\t0\tflip_swap"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert read_bed_matrix(out_prefix.with_suffix(".bed"), 2, 1) == [[3, 0]]


def test_apply_vmap_to_bfile_reports_chunk_progress(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(
        source_prefix,
        ["1\trs1\t0\t100\tA\tG", "1\trs2\t0\t200\tC\tT"],
        ["S1 S1 0 0 0 -9"],
        [[0], [3]],
    )
    write_lines(
        vmap,
        [
            "1\t100\tt1\tA\tG\t.\t0\tidentity",
            "1\t200\tt2\tC\tT\t.\t1\tidentity",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py_with_env(
        "apply_vmap_to_bfile.py",
        {"MATCH_BFILE_APPLY_CHUNK_SIZE": "1"},
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
    )
    assert result.returncode == 0, result.stderr
    assert "apply_vmap_to_bfile.py progress: 1/2 chunks" in result.stderr
    assert "apply_vmap_to_bfile.py progress: 2/2 chunks" in result.stderr


def test_apply_vmap_to_bfile_chrX_chunk_boundary_parity(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_chunk1 = tmp_path / "aligned_chunk1"
    out_chunk3 = tmp_path / "aligned_chunk3"
    # Position 5,000,000 on chrX is outside PAR for GRCh37, so male ploidy differs from female.
    write_bfile(
        source_prefix,
        [
            "X\trs1\t0\t5000000\tA\tG",
            "X\trs2\t0\t5000001\tC\tT",
            "X\trs3\t0\t5000002\tG\tA",
        ],
        ["S1 S1 0 0 1 -9", "S2 S2 0 0 2 -9", "S3 S3 0 0 0 -9"],
        [[2, 3, 0], [0, 0, 3], [3, 2, 1]],
    )
    write_lines(
        vmap,
        [
            "X\t5000000\tt1\tA\tG\t.\t0\tidentity",
            "X\t5000001\tt2\tC\tT\t.\t1\tidentity",
            "X\t5000002\tt3\tG\tA\t.\t2\tidentity",
        ],
    )
    write_json(
        vmap.with_name(vmap.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result_chunk1 = run_py_with_env(
        "apply_vmap_to_bfile.py",
        {"MATCH_BFILE_APPLY_CHUNK_SIZE": "1"},
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_chunk1,
    )
    result_chunk3 = run_py_with_env(
        "apply_vmap_to_bfile.py",
        {"MATCH_BFILE_APPLY_CHUNK_SIZE": "3"},
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_chunk3,
    )

    assert result_chunk1.returncode == 0, result_chunk1.stderr
    assert result_chunk3.returncode == 0, result_chunk3.stderr

    assert read_bed_matrix(out_chunk1.with_suffix(".bed"), 3, 3) == read_bed_matrix(out_chunk3.with_suffix(".bed"), 3, 3)
    assert out_chunk1.with_suffix(".bim").read_text(encoding="utf-8") == out_chunk3.with_suffix(".bim").read_text(encoding="utf-8")
    assert out_chunk1.with_suffix(".fam").read_text(encoding="utf-8") == out_chunk3.with_suffix(".fam").read_text(encoding="utf-8")

    expected_haploid_warning = "found 1 incompatible heterozygous haploid-target genotypes"
    expected_unknown_warning = "skipped ploidy validation for 3 sex-dependent row/sample cells with unknown output sex"
    assert expected_haploid_warning in result_chunk1.stderr
    assert expected_haploid_warning in result_chunk3.stderr
    assert expected_unknown_warning in result_chunk1.stderr
    assert expected_unknown_warning in result_chunk3.stderr


def test_apply_vmap_to_bfile_warns_about_unknown_sex_only_for_sex_dependent_rows(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["X\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 0 -9", "S2 S2 0 0 1 -9"], [[0, 3]])
    write_lines(vmap, ["X\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert "skipped ploidy validation for 1 sex-dependent row/sample cells with unknown output sex" in result.stderr
    assert read_bed_matrix(out_prefix.with_suffix(".bed"), 2, 1) == [[0, 3]]
    assert out_prefix.with_suffix(".ploidy").exists()

def test_apply_vmap_to_bfile_validates_sex_independent_haploid_rows_even_with_unknown_sex(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["MT\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 0 -9"], [[2]])
    write_lines(vmap, ["MT\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert "found 1 incompatible heterozygous haploid-target genotypes" in result.stderr
    assert read_bed_matrix(out_prefix.with_suffix(".bed"), 1, 1) == [[2]]
    assert out_prefix.with_suffix(".ploidy").exists()


def test_apply_vmap_to_bfile_reports_haploid_heterozygous_calls_without_rewriting(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["X\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 1 -9"], [[2]])
    write_lines(vmap, ["X\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert "found 1 incompatible heterozygous haploid-target genotypes" in result.stderr
    assert read_bed_matrix(out_prefix.with_suffix(".bed"), 1, 1) == [[2]]
    assert out_prefix.with_suffix(".ploidy").exists()


def test_apply_vmap_to_bfile_reports_absent_region_calls_without_rewriting(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["Y\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 2 -9"], [[0]])
    write_lines(vmap, ["Y\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert "found 1 incompatible absent-region nonmissing genotypes" in result.stderr
    assert read_bed_matrix(out_prefix.with_suffix(".bed"), 1, 1) == [[0]]
    assert out_prefix.with_suffix(".ploidy").exists()


def test_apply_vmap_to_bfile_writes_qc_tsv_for_haploid_het_incompatible(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    # Male with het genotype at X non-PAR (X:100 is outside PAR1 for GRCh37) → haploid_het_incompatible
    write_bfile(source_prefix, ["X\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 1 -9"], [[2]])
    write_lines(vmap, ["X\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr

    qc_path = out_prefix.with_suffix(".qc.tsv")
    assert qc_path.exists()
    lines = qc_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "source_shard\tsource_index\tid\tstatus\tn_affected"
    assert lines[1] == ".\t0\tX:100:A:G\thaploid_het_incompatible\t1"
    assert len(lines) == 2


def test_apply_vmap_to_bfile_writes_qc_tsv_for_absent_nonmissing(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    # Female with non-missing genotype at Y MSY (Y:100 is in Y_MSY for GRCh37) → absent_nonmissing
    write_bfile(source_prefix, ["Y\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 2 -9"], [[0]])
    write_lines(vmap, ["Y\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr

    qc_path = out_prefix.with_suffix(".qc.tsv")
    assert qc_path.exists()
    lines = qc_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "source_shard\tsource_index\tid\tstatus\tn_affected"
    assert lines[1] == ".\t0\tY:100:A:G\tabsent_nonmissing\t1"
    assert len(lines) == 2


def test_apply_vmap_to_bfile_no_qc_tsv_when_no_ploidy_issues(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 1 -9"], [[0]])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert not out_prefix.with_suffix(".qc.tsv").exists()


def test_apply_vmap_to_bfile_removes_stale_qc_tsv_when_no_ploidy_issues(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 1 -9"], [[0]])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})
    out_prefix.with_suffix(".qc.tsv").write_text("stale content\n", encoding="utf-8")

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert not out_prefix.with_suffix(".qc.tsv").exists()


def test_apply_vmap_to_bfile_uses_exact_source_shard_for_sharded_payloads(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    fam_lines = ["S1 S1 0 0 0 -9"]
    write_bfile(tmp_path / "source_chr1", ["chr1\trs1\t0\t100\tA\tG"], fam_lines, [[0]])
    write_bfile(tmp_path / "source_2", ["2\trs2\t0\t200\tC\tT"], fam_lines, [[3]])
    for stem in ("source_chr1", "source_2"):
        src = tmp_path / stem
        token = stem.split("_", 1)[1]
        for ext in (".bed", ".bim", ".fam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_lines(vmap, ["1\t100\tt1\tA\tG\tchr1\t0\tidentity", "2\t200\tt2\tC\tT\t2\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_template, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert [row.snp for row in read_bim(out_prefix.with_suffix(".bim"))] == ["1:100:A:G", "2:200:C:T"]


def test_apply_vmap_to_bfile_discovers_xy_shard_alias_for_chr_template(tmp_path):
    source_template = tmp_path / "source_chr@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    fam_lines = ["S1 S1 0 0 2 -9"]
    write_bfile(tmp_path / "source_chrX", ["23\trs_x\t0\t100\tA\tG"], fam_lines, [[0]])
    write_bfile(tmp_path / "source_chrXY", ["25\trs_xy\t0\t200\tC\tT"], fam_lines, [[3]])
    write_lines(vmap, ["X\t100\tt_x\tA\tG\tX\t0\tidentity", "X\t200\tt_xy\tC\tT\tXY\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_template, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert [row.snp for row in read_bim(out_prefix.with_suffix(".bim"))] == ["X:100:A:G", "X:200:C:T"]
    assert read_bed_matrix(out_prefix.with_suffix(".bed"), 1, 2) == [[0], [3]]


def test_apply_vmap_to_bfile_rejects_missing_required_source_shard(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(tmp_path / "source_chr1", ["chr1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 0 -9"], [[0]])
    for ext in (".bed", ".bim", ".fam"):
        (tmp_path / f"source_chr1{ext}").rename(tmp_path / f"source.chr1{ext}")
    write_lines(vmap, ["1\t100\tt1\tA\tG\t1\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_template, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode != 0
    assert "missing required source shard" in result.stderr


def test_apply_vmap_to_bfile_rejects_zero_discovered_source_shards(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_lines(vmap, ["1\t100\tt1\tA\tG\t1\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_template, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode != 0
    assert "no source shards found for template" in result.stderr


def test_apply_vmap_to_bfile_rejects_all_unmatched_target_rows(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 0 -9"], [[0]])
    write_lines(vmap, ["1\t400\tt_missing\tA\tC\t.\t-1\tmissing"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode != 0
    assert "only unmatched target rows" in result.stderr


def test_apply_vmap_to_bfile_rejects_empty_output_after_only_mapped_target_filter(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 0 -9"], [[0]])
    write_lines(vmap, ["1\t400\tt_missing\tA\tC\t.\t-1\tmissing"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py(
        "apply_vmap_to_bfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--only-mapped-target",
    )

    assert result.returncode != 0
    assert "no retained target rows remain" in result.stderr


def test_apply_vmap_to_bfile_shards_output_by_target_contig(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned.@"
    write_bfile(source_prefix, ["1\trs1\t0\t100\tA\tG", "2\trs2\t0\t200\tC\tT"], ["S1 S1 0 0 0 -9"], [[0], [3]])
    write_lines(vmap, ["2\t200\tt2\tC\tT\t.\t1\tidentity", "1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "aligned.2.bim").exists()
    assert (tmp_path / "aligned.1.bim").exists()
    assert (tmp_path / "aligned.2.fam").read_text(encoding="utf-8") == source_prefix.with_suffix(".fam").read_text(encoding="utf-8")
    assert (tmp_path / "aligned.1.fam").read_text(encoding="utf-8") == source_prefix.with_suffix(".fam").read_text(encoding="utf-8")


def test_apply_vmap_to_bfile_suppresses_empty_output_shards(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned.@"
    write_bfile(source_prefix, ["2\trs2\t0\t200\tC\tT"], ["S1 S1 0 0 0 -9"], [[3]])
    write_lines(vmap, ["2\t200\tt2\tC\tT\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "aligned.2.bim").exists()
    assert not (tmp_path / "aligned.1.bim").exists()


def test_apply_vmap_to_bfile_fails_on_out_of_range_single_file_provenance(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 0 -9"], [[0]])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t5\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode != 0
    assert "source provenance out of range" in result.stderr


def test_apply_vmap_to_bfile_supports_dotted_output_prefix(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned.panel"
    write_bfile(source_prefix, ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 0 -9"], [[0]])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    bim_path = out_prefix.parent / f"{out_prefix.name}.bim"
    fam_path = out_prefix.parent / f"{out_prefix.name}.fam"
    assert [(row.chrom, row.snp, row.cm, row.bp, row.a1, row.a2) for row in read_bim(bim_path)] == [
        ("1", "1:100:A:G", "0", "100", "A", "G"),
    ]
    assert read_fam_sample_count(fam_path) == 1


def test_apply_vmap_to_bfile_merges_source_shards_when_target_contigs_merge(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned.@"
    fam_lines = ["S1 S1 0 0 2 -9", "S2 S2 0 0 2 -9"]
    write_bfile(tmp_path / "source_23", ["23\trs_nonpar\t0\t70000\tA\tG"], fam_lines, [[0, 3]])
    write_bfile(tmp_path / "source_25", ["25\trs_par\t0\t60001\tC\tT"], fam_lines, [[2, 0]])
    for stem in ("source_23", "source_25"):
        src = tmp_path / stem
        token = stem.split("_", 1)[1]
        for ext in (".bed", ".bim", ".fam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_lines(
        vmap,
        [
            "X\t70000\tt_nonpar\tA\tG\t23\t0\tidentity",
            "X\t60001\tt_par\tC\tT\t25\t0\tidentity",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_template, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert [(row.chrom, row.snp) for row in read_bim(tmp_path / "aligned.X.bim")] == [("X", "X:70000:A:G"), ("X", "X:60001:C:T")]
    assert read_bed_matrix(tmp_path / "aligned.X.bed", 2, 2) == [[0, 3], [2, 0]]


def test_apply_vmap_to_bfile_splits_source_shard_when_target_contigs_split(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned.@"
    fam_lines = ["S1 S1 0 0 2 -9", "S2 S2 0 0 2 -9"]
    write_bfile(tmp_path / "source_23", ["23\trs_nonpar\t0\t70000\tA\tG", "23\trs_par\t0\t60001\tC\tT"], fam_lines, [[0, 3], [3, 0]])
    for ext in (".bed", ".bim", ".fam"):
        (tmp_path / f"source_23{ext}").rename(tmp_path / f"source.23{ext}")
    write_lines(
        vmap,
        [
            "25\t60001\tt_par\tC\tT\t23\t1\tidentity",
            "23\t70000\tt_nonpar\tA\tG\t23\t0\tidentity",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "plink_splitx"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_template, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert [(row.chrom, row.snp) for row in read_bim(tmp_path / "aligned.25.bim")] == [("25", "25:60001:C:T")]
    assert [(row.chrom, row.snp) for row in read_bim(tmp_path / "aligned.23.bim")] == [("23", "23:70000:A:G")]
    assert read_bed_matrix(tmp_path / "aligned.25.bed", 2, 1) == [[3, 0]]
    assert read_bed_matrix(tmp_path / "aligned.23.bed", 2, 1) == [[0, 3]]


def test_apply_vmap_to_bfile_target_fam_reorders_and_fills_missing_for_inconsistent_shards(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    target_fam = tmp_path / "target.fam"
    write_bfile(tmp_path / "source_1", ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 1 -9", "S2 S2 0 0 2 -9"], [[0, 3]])
    write_bfile(tmp_path / "source_2", ["2\trs2\t0\t200\tC\tT"], ["S2 S2 0 0 2 -9", "S3 S3 0 0 1 -9"], [[2, 0]])
    for stem in ("source_1", "source_2"):
        src = tmp_path / stem
        token = stem.split("_", 1)[1]
        for ext in (".bed", ".bim", ".fam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_lines(vmap, ["1\t100\tt1\tA\tG\t1\t0\tidentity", "2\t200\tt2\tC\tT\t2\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})
    write_lines(target_fam, ["S3 S3 0 0 1 -9", "S1 S1 0 0 1 -9", "S2 S2 0 0 2 -9"])

    result = run_py(
        "apply_vmap_to_bfile.py",
        "--source-prefix",
        source_template,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--target-fam",
        target_fam,
    )

    assert result.returncode == 0, result.stderr
    assert out_prefix.with_suffix(".fam").read_text(encoding="utf-8") == target_fam.read_text(encoding="utf-8")
    assert read_bed_matrix(out_prefix.with_suffix(".bed"), 3, 2) == [[1, 0, 3], [0, 1, 2]]
    assert "Sample-axis reconciliation summary" in result.stderr


def test_apply_vmap_to_bfile_iid_mode_ignores_fid_for_target_fam(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    target_fam = tmp_path / "target.fam"
    write_bfile(source_prefix, ["1\trs1\t0\t100\tA\tG"], ["F1 S1 0 0 1 -9", "F2 S2 0 0 2 -9"], [[0, 3]])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})
    write_lines(target_fam, ["X9 S2 0 0 2 -9", "X8 S1 0 0 1 -9"])

    result = run_py(
        "apply_vmap_to_bfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--target-fam",
        target_fam,
        "--sample-id-mode",
        "iid",
    )

    assert result.returncode == 0, result.stderr
    assert read_bed_matrix(out_prefix.with_suffix(".bed"), 2, 1) == [[3, 0]]


def test_apply_vmap_to_bfile_rejects_duplicate_subject_keys_in_source_shard(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["1\trs1\t0\t100\tA\tG"], ["F1 S1 0 0 1 -9", "F2 S1 0 0 2 -9"], [[0, 3]])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py(
        "apply_vmap_to_bfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--sample-id-mode",
        "iid",
    )

    assert result.returncode != 0
    assert "duplicate subject key" in result.stderr


def test_apply_vmap_to_bfile_target_fam_controls_haploid_validation_sex(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    target_fam = tmp_path / "target.fam"
    write_bfile(source_prefix, ["X\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 0 -9"], [[0]])
    write_lines(vmap, ["X\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})
    write_lines(target_fam, ["S1 S1 0 0 1 -9"])

    result = run_py(
        "apply_vmap_to_bfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--target-fam",
        target_fam,
    )

    assert result.returncode == 0, result.stderr
    assert "skipping haploid validation" not in result.stderr
    assert out_prefix.with_suffix(".ploidy").exists()


def test_apply_vmap_to_bfile_relocates_genotypes_across_output_shards_by_provenance(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned.@"
    fam_lines = ["S1 S1 0 0 0 -9", "S2 S2 0 0 0 -9"]
    write_bfile(tmp_path / "source_1", ["1\trs1\t0\t100\tA\tG"], fam_lines, [[0, 3]])
    write_bfile(tmp_path / "source_2", ["2\trs2\t0\t200\tC\tT"], fam_lines, [[3, 0]])
    for stem in ("source_1", "source_2"):
        src = tmp_path / stem
        token = stem.split("_", 1)[1]
        for ext in (".bed", ".bim", ".fam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_lines(
        vmap,
        [
            "2\t500\tt_from_1\tA\tG\t1\t0\tidentity",
            "1\t600\tt_from_2\tC\tT\t2\t0\tidentity",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_template, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode == 0, result.stderr
    assert [(row.chrom, row.snp) for row in read_bim(tmp_path / "aligned.2.bim")] == [("2", "2:500:A:G")]
    assert [(row.chrom, row.snp) for row in read_bim(tmp_path / "aligned.1.bim")] == [("1", "1:600:C:T")]
    assert read_bed_matrix(tmp_path / "aligned.2.bed", 2, 1) == [[0, 3]]
    assert read_bed_matrix(tmp_path / "aligned.1.bed", 2, 1) == [[3, 0]]
    assert (tmp_path / "aligned.2.fam").read_text(encoding="utf-8") == (tmp_path / "source.1.fam").read_text(encoding="utf-8")
    assert (tmp_path / "aligned.1.fam").read_text(encoding="utf-8") == (tmp_path / "source.1.fam").read_text(encoding="utf-8")


def test_apply_vmap_to_bfile_rejects_sharded_source_fam_content_mismatch(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(tmp_path / "source_1", ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 1 -9"], [[0]])
    write_bfile(tmp_path / "source_2", ["2\trs2\t0\t200\tC\tT"], ["S1 S1 0 0 1 1"], [[3]])
    for stem in ("source_1", "source_2"):
        src = tmp_path / stem
        token = stem.split("_", 1)[1]
        for ext in (".bed", ".bim", ".fam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_lines(vmap, ["1\t100\tt1\tA\tG\t1\t0\tidentity", "2\t200\tt2\tC\tT\t2\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_bfile.py", "--source-prefix", source_template, "--vmap", vmap, "--output-prefix", out_prefix)
    assert result.returncode != 0
    assert "identical .fam contents" in result.stderr


def test_apply_vmap_to_bfile_sample_axis_native_keeps_per_output_shard_fam(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned.@"
    write_bfile(tmp_path / "source_1", ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 1 -9", "S2 S2 0 0 2 -9"], [[0, 3]])
    write_bfile(tmp_path / "source_X", ["X\trsx\t0\t5000000\tC\tT"], ["S1 S1 0 0 1 -9"], [[2]])
    for stem, token in (("source_1", "1"), ("source_X", "X")):
        src = tmp_path / stem
        for ext in (".bed", ".bim", ".fam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_lines(vmap, ["1\t100\tt1\tA\tG\t1\t0\tidentity", "X\t5000000\ttx\tC\tT\tX\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py(
        "apply_vmap_to_bfile.py",
        "--source-prefix",
        source_template,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--sample-axis",
        "native",
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "aligned.1.fam").read_text(encoding="utf-8") == (tmp_path / "source.1.fam").read_text(encoding="utf-8")
    assert (tmp_path / "aligned.X.fam").read_text(encoding="utf-8") == (tmp_path / "source.X.fam").read_text(encoding="utf-8")
    assert read_bed_matrix(tmp_path / "aligned.1.bed", 2, 1) == [[0, 3]]
    assert read_bed_matrix(tmp_path / "aligned.X.bed", 1, 1) == [[2]]


def test_apply_vmap_to_bfile_sample_axis_native_rejects_mixed_axes_in_one_output(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(tmp_path / "source_1", ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 1 -9", "S2 S2 0 0 2 -9"], [[0, 3]])
    write_bfile(tmp_path / "source_X", ["X\trsx\t0\t5000000\tC\tT"], ["S1 S1 0 0 1 -9"], [[2]])
    for stem, token in (("source_1", "1"), ("source_X", "X")):
        src = tmp_path / stem
        for ext in (".bed", ".bim", ".fam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_lines(vmap, ["1\t100\tt1\tA\tG\t1\t0\tidentity", "X\t5000000\ttx\tC\tT\tX\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py(
        "apply_vmap_to_bfile.py",
        "--source-prefix",
        source_template,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--sample-axis",
        "native",
    )

    assert result.returncode != 0
    assert "--sample-axis native cannot emit" in result.stderr
    assert "different .fam contents" in result.stderr


def test_apply_vmap_to_bfile_skip_ploidy_check_suppresses_validation_warning_and_qc(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(source_prefix, ["X\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 1 -9"], [[2]])
    write_lines(vmap, ["X\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py(
        "apply_vmap_to_bfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--skip-ploidy-check",
    )

    assert result.returncode == 0, result.stderr
    assert "incompatible heterozygous haploid-target genotypes" not in result.stderr
    assert not out_prefix.with_suffix(".qc.tsv").exists()
    assert out_prefix.with_suffix(".ploidy").exists()


def test_apply_vmap_to_bfile_streams_single_file_payload_in_small_chunks(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_bfile(
        source_prefix,
        ["1\trs1\t0\t100\tA\tG", "1\trs2\t0\t200\tC\tT", "1\trs3\t0\t300\tG\tA"],
        ["S1 S1 0 0 0 -9", "S2 S2 0 0 0 -9"],
        [[0, 3], [0, 2], [3, 0]],
    )
    write_lines(
        vmap,
        [
            "1\t200\tt2\tT\tC\t.\t1\tswap",
            "1\t100\tt1\tA\tG\t.\t0\tidentity",
            "1\t400\tt3\tA\tC\t.\t-1\tmissing",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py_with_env(
        "apply_vmap_to_bfile.py",
        {"MATCH_BFILE_APPLY_CHUNK_SIZE": "1"},
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
    )

    assert result.returncode == 0, result.stderr
    assert [(row.chrom, row.snp) for row in read_bim(out_prefix.with_suffix(".bim"))] == [("1", "1:200:T:C"), ("1", "1:100:A:G"), ("1", "1:400:A:C")]
    assert read_bed_matrix(out_prefix.with_suffix(".bed"), 2, 3) == [[3, 2], [0, 3], [1, 1]]


def test_apply_vmap_to_bfile_streams_interleaved_sharded_input_with_only_mapped_target(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    fam_lines = ["S1 S1 0 0 0 -9", "S2 S2 0 0 0 -9"]
    write_bfile(tmp_path / "source_1", ["1\trs1\t0\t100\tA\tG", "1\trs1b\t0\t150\tC\tT"], fam_lines, [[0, 3], [2, 0]])
    write_bfile(tmp_path / "source_2", ["2\trs2\t0\t200\tC\tT"], fam_lines, [[3, 0]])
    for stem in ("source_1", "source_2"):
        src = tmp_path / stem
        token = stem.split("_", 1)[1]
        for ext in (".bed", ".bim", ".fam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_lines(
        vmap,
        [
            "2\t500\tt_from_1a\tA\tG\t1\t0\tidentity",
            "1\t600\tt_from_2\tC\tT\t2\t0\tidentity",
            "1\t700\tt_missing\tA\tC\t.\t-1\tmissing",
            "2\t800\tt_from_1b\tT\tC\t1\t1\tswap",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py_with_env(
        "apply_vmap_to_bfile.py",
        {"MATCH_BFILE_APPLY_CHUNK_SIZE": "1"},
        "--source-prefix",
        source_template,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--only-mapped-target",
    )

    assert result.returncode == 0, result.stderr
    assert [(row.chrom, row.snp) for row in read_bim(out_prefix.with_suffix(".bim"))] == [
        ("2", "2:500:A:G"),
        ("1", "1:600:C:T"),
        ("2", "2:800:T:C"),
    ]
    assert read_bed_matrix(out_prefix.with_suffix(".bed"), 2, 3) == [[0, 3], [3, 0], [2, 3]]

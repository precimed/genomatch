import gzip
import math

from utils import run_py, write_json, write_lines


def write_vmap_with_meta(path, lines):
    write_lines(path, lines)
    write_json(
        path.with_name(path.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )


def test_apply_vmap_to_sumstats_target_order_and_swap_effects(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA\tOR\tEAF", "1\t100\tA\tG\t0.5\t2.0\t0.1", "1\t200\tC\tT\t1.0\t4.0\t0.3"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA", "col_OR: OR", "col_EAF: EAF"])
    write_lines(
        vmap,
        [
            "1\t200\tt2\tT\tC\t.\t1\tswap",
            "1\t100\tt1\tA\tG\t.\t0\tidentity",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tEA\tOA\tBETA\tOR\tEAF\n"
        "1\t200\tT\tC\t-1.0\t0.25\t0.7\n"
        "1\t100\tA\tG\t0.5\t2.0\t0.1\n"
    )


def test_apply_vmap_to_sumstats_appends_missing_variant_columns(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tBETA", "1\t0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_SNP: SNP", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR\tBETA\tPOS\tSNP\tEA\tOA\n1\t0.5\t100\tt1\tA\tG\n"


def test_apply_vmap_to_sumstats_keeps_unmatched_target_rows_by_default(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA\tG\t0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity", "1\t200\tt2\tC\tT\t.\t-1\tmissing"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tEA\tOA\tBETA\n1\t100\tA\tG\t0.5\n1\t200\tC\tT\tn/a\n"


def test_apply_vmap_to_sumstats_only_mapped_target_drops_unmatched_rows(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA\tG\t0.5", "1\t300\tG\tA\t1.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(
        vmap,
        [
            "1\t200\tt_missing\tC\tT\t.\t-1\tmissing",
            "1\t300\tt3\tG\tA\t.\t1\tidentity",
            "1\t100\tt1\tA\tG\t.\t0\tidentity",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out,
        "--only-mapped-target",
    )
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tEA\tOA\tBETA\n1\t300\tG\tA\t1.5\n1\t100\tA\tG\t0.5\n"


def test_apply_vmap_to_sumstats_warns_and_writes_na_for_invalid_swapped_numeric_effects(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA\tOR\tEAF\tORL95\tORU95", "1\t100\tA\tG\tbad\t0\tnope\tnan\tinf"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_BETA: BETA",
            "col_OR: OR",
            "col_EAF: EAF",
            "col_ORL95: ORL95",
            "col_ORU95: ORU95",
        ],
    )
    write_lines(vmap, ["1\t100\tt1\tG\tA\t.\t0\tswap"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert "writing n/a" in result.stderr.lower()
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tEA\tOA\tBETA\tOR\tEAF\tORL95\tORU95\n1\t100\tG\tA\tn/a\tn/a\tn/a\tn/a\tn/a\n"


def test_apply_vmap_to_sumstats_flip_swap_uses_swap_numeric_effects(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA\tOR\tEAF", "1\t100\tA\tC\t0.5\t2.0\t0.1"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA", "col_OR: OR", "col_EAF: EAF"])
    write_lines(vmap, ["1\t100\tt1\tG\tT\t.\t0\tflip_swap"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tEA\tOA\tBETA\tOR\tEAF\n1\t100\tG\tT\t-0.5\t0.5\t0.9\n"


def test_apply_vmap_to_sumstats_rewrites_target_chr_pos_and_id_when_columns_exist(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t100\trs_source\tA\tG\t0.5"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_BETA: BETA",
        ],
    )
    write_lines(vmap, ["2\t250\trs_target\tT\tC\t.\t0\tswap"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEA\tOA\tBETA\n2\t250\trs_target\tT\tC\t-0.5\n"


def test_apply_vmap_to_sumstats_swap_inverts_and_swaps_or_confidence_bounds(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tORL95\tORU95", "1\t100\tA\tG\t1.2\t1.5"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_ORL95: ORL95",
            "col_ORU95: ORU95",
        ],
    )
    write_lines(vmap, ["1\t100\tt1\tG\tA\t.\t0\tswap"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == f"CHR\tPOS\tEA\tOA\tORL95\tORU95\n1\t100\tG\tA\t{1.0 / 1.5}\t{1.0 / 1.2}\n"


def test_apply_vmap_to_sumstats_reads_joined_source_variant_fields(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["VARIANT\tBETA", "1_100:A:G\t0.5"])
    write_lines(
        meta,
        [
            "col_CHR: VARIANT",
            "col_POS: VARIANT",
            "col_EffectAllele: VARIANT",
            "col_OtherAllele: VARIANT",
            "col_BETA: BETA",
        ],
    )
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "VARIANT\tBETA\n1_100:A:G\t0.5\n"


def test_apply_vmap_to_sumstats_rewrites_joined_chr_field_in_place(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["VARIANT\tBETA", "1_100:A:G\t0.5"])
    write_lines(meta, ["col_CHR: VARIANT", "col_BETA: BETA"])
    write_lines(vmap, ["2\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "VARIANT\tBETA\n2_100:A:G\t0.5\n"


def test_apply_vmap_to_sumstats_preserves_csv_delimiter_in_output(tmp_path):
    sumstats = tmp_path / "ss.csv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.csv"
    write_lines(sumstats, ["CHR,POS,EA,OA,BETA", "1,100,A,G,0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(vmap, ["1\t100\tt1\tG\tA\t.\t0\tswap"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR,POS,EA,OA,BETA\n1,100,G,A,-0.5\n"


def test_apply_vmap_to_sumstats_writes_gzipped_output_when_suffix_is_gz(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv.gz"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA\tG\t0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(vmap, ["1\t100\tt1\tG\tA\t.\t0\tswap"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    with gzip.open(out, "rt", encoding="utf-8") as handle:
        assert handle.read() == "CHR\tPOS\tEA\tOA\tBETA\n1\t100\tG\tA\t-0.5\n"


def test_apply_vmap_to_sumstats_fails_on_out_of_range_single_file_provenance(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA\tG\t0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t5\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode != 0
    assert "source_shard='.'" in result.stderr


def test_apply_vmap_to_sumstats_rejects_non_single_file_provenance(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA\tG\t0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t1\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode != 0
    assert "single-file payload lookup only for source_shard='.'" in result.stderr


def test_apply_vmap_to_sumstats_rejects_template_input_path(tmp_path):
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", tmp_path / "ss.@.tsv", "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode != 0
    assert "does not accept '@' paths" in result.stderr


def test_apply_vmap_to_sumstats_rejects_template_output_path(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA\tG\t0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", tmp_path / "out.@.tsv")
    assert result.returncode != 0
    assert "does not accept '@' paths" in result.stderr


def test_apply_vmap_to_sumstats_fails_on_short_payload_rows(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode != 0
    assert "fewer columns than expected" in result.stderr


def test_apply_vmap_to_sumstats_clean_normalizes_headers_and_drops_unrecognized_columns(tmp_path):
    sumstats = tmp_path / "ss.csv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(
        sumstats,
        [
            "CHR,POS,SNP,EA,OA,p-value,Beta!,unused,INFO??",
            "1,100,rs1,A,G,0.05,0.5,dropme,",
        ],
    )
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_P: p value",
            "col_BETA: beta!",
            "col_INFO: info??",
            "stats_Model: linear",
        ],
    )
    write_vmap_with_meta(vmap, ["1\t100\trs_target\tT\tC\t.\t0\tidentity"])

    result = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out,
        "--clean",
    )

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\tZ\tBETA\tSE\n"
        "1\t100\trs_target\tT\tC\t0.05\t1.95996398454005\t0.5\t0.255106728462327\n"
    )


def test_apply_vmap_to_sumstats_clean_fails_on_ambiguous_normalized_header_lookup(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tSNP\tEA\tOA\tP!\tP", "1\t100\trs1\tA\tG\t0.05\t0.06"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_P: P",
            "stats_Model: linear",
        ],
    )
    write_vmap_with_meta(vmap, ["1\t100\trs1\tA\tG\t.\t0\tidentity"])

    result = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out,
        "--clean",
    )

    assert result.returncode != 0
    assert "identified uniquely" in result.stderr


def test_apply_vmap_to_sumstats_clean_converts_log_p_variants(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out_neg = tmp_path / "out-neg.tsv"
    out_log = tmp_path / "out-log.tsv"
    write_lines(sumstats, ["CHR\tPOS\tSNP\tEA\tOA\tPV", "1\t100\trs1\tA\tG\t2"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_P: PV",
            "stats_Model: linear",
            "stats_neglog10P: true",
        ],
    )
    write_vmap_with_meta(vmap, ["1\t100\trs1\tA\tG\t.\t0\tidentity"])

    neg = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out_neg,
        "--clean",
    )
    assert neg.returncode == 0, neg.stderr
    assert out_neg.read_text(encoding="utf-8").endswith("\t0.01\n")

    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_P: PV",
            "stats_Model: linear",
            "stats_log10P: true",
        ],
    )
    log_result = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out_log,
        "--clean",
    )
    assert log_result.returncode == 0, log_result.stderr
    assert out_log.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\n1\t100\trs1\tA\tG\t\n"


def test_apply_vmap_to_sumstats_clean_fill_mode_column_and_row_differ(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out_column = tmp_path / "out-column.tsv"
    out_row = tmp_path / "out-row.tsv"
    write_lines(
        sumstats,
        [
            "CHR\tPOS\tSNP\tEA\tOA\tBETA\tSE\tZ",
            "1\t100\trs1\tA\tG\t2\t1\t",
            "1\t200\trs2\tC\tT\t4\t2\t9",
        ],
    )
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_BETA: BETA",
            "col_SE: SE",
            "col_Z: Z",
            "stats_Model: linear",
        ],
    )
    write_vmap_with_meta(vmap, ["1\t100\trs1\tA\tG\t.\t0\tidentity", "1\t200\trs2\tC\tT\t.\t1\tidentity"])

    column = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out_column,
        "--clean",
    )
    assert column.returncode == 0, column.stderr
    assert out_column.read_text(encoding="utf-8").splitlines() == [
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\tZ\tBETA\tSE",
        "1\t100\trs1\tA\tG\t\t\t2\t1",
        "1\t200\trs2\tC\tT\t0\t9\t4\t2",
    ]

    row = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out_row,
        "--clean",
        "--fill-mode",
        "row",
    )
    assert row.returncode == 0, row.stderr
    assert out_row.read_text(encoding="utf-8").splitlines() == [
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\tZ\tBETA\tSE",
        "1\t100\trs1\tA\tG\t0.0455002638963584\t2\t2\t1",
        "1\t200\trs2\tC\tT\t0\t9\t4\t2",
    ]


def test_apply_vmap_to_sumstats_clean_derives_eaf_from_oaf_and_drops_oaf_columns(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tSNP\tEA\tOA\tOAF", "1\t100\trs1\tA\tG\t0.2"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_OAF: OAF",
            "stats_Model: linear",
        ],
    )
    write_vmap_with_meta(vmap, ["1\t100\trs1\tA\tG\t.\t0\tidentity"])

    result = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out,
        "--clean",
    )

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tEAF\n"
        "1\t100\trs1\tA\tG\t0.8\n"
    )


def test_apply_vmap_to_sumstats_clean_prefers_z_from_p_and_beta_and_runs_second_validation(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(
        sumstats,
        [
            "CHR\tPOS\tSNP\tEA\tOA\tP\tBETA\tSE\tORL95\tORU95",
            "1\t100\trs1\tA\tG\t0.0455002638963584\t1\t1\t2\t1",
        ],
    )
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_P: P",
            "col_BETA: BETA",
            "col_SE: SE",
            "col_ORL95: ORL95",
            "col_ORU95: ORU95",
            "stats_Model: logistic",
        ],
    )
    write_vmap_with_meta(vmap, ["1\t100\trs1\tA\tG\t.\t0\tidentity"])

    result = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out,
        "--clean",
        "--fill-mode",
        "row",
    )

    assert result.returncode == 0, result.stderr
    lines = out.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    row = lines[1].split("\t")
    values = dict(zip(header, row))
    assert math.isclose(float(values["Z"]), 2.0, rel_tol=1e-9)
    assert values["SE"] == "1"


def test_apply_vmap_to_sumstats_clean_use_af_inference_derives_n(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tSNP\tEA\tOA\tZ\tBETA\tEAF", "1\t100\trs1\tA\tG\t2\t0.5\t0.25"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_Z: Z",
            "col_BETA: BETA",
            "col_EAF: EAF",
            "stats_Model: linear",
        ],
    )
    write_vmap_with_meta(vmap, ["1\t100\trs1\tA\tG\t.\t0\tidentity"])

    result = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out,
        "--clean",
        "--use-af-inference",
    )

    assert result.returncode == 0, result.stderr
    header, row = out.read_text(encoding="utf-8").splitlines()
    values = dict(zip(header.split("\t"), row.split("\t")))
    assert math.isclose(float(values["N"]), 38.6666666666667, rel_tol=1e-12)


def test_apply_vmap_to_sumstats_clean_keeps_unmatched_rows_all_missing(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(
        sumstats,
        [
            "CHR\tPOS\tSNP\tEA\tOA\tBETA",
            "1\t100\trs1\tA\tG\t0.5",
        ],
    )
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_BETA: BETA",
            "stats_Model: linear",
            "stats_TotalN: 100",
        ],
    )
    write_vmap_with_meta(vmap, ["1\t100\trs1\tA\tG\t.\t0\tidentity", "1\t200\trs2\tC\tT\t.\t-1\tmissing"])

    result = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out,
        "--clean",
    )

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tN\tBETA\n"
        "1\t100\trs1\tA\tG\t100\t0.5\n"
        "1\t200\trs2\tC\tT\t\t\n"
    )


def test_apply_vmap_to_sumstats_clean_swap_keeps_direction_and_drops_missing_p_under_only_mapped_target(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(
        sumstats,
        [
            "CHR\tPOS\tSNP\tEA\tOA\tP\tBETA\tEAF\tDirection",
            "1\t100\trs1\tA\tG\t0.05\t0.5\t0.2\t+-",
            "1\t200\trs2\tC\tT\t\t1.5\t0.4\t--",
        ],
    )
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_P: P",
            "col_BETA: BETA",
            "col_EAF: EAF",
            "col_Direction: Direction",
            "stats_Model: linear",
        ],
    )
    write_vmap_with_meta(vmap, ["1\t100\trs1\tG\tA\t.\t0\tswap", "1\t200\trs2\tT\tC\t.\t1\tidentity"])

    result = run_py(
        "apply_vmap_to_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--vmap",
        vmap,
        "--output",
        out,
        "--clean",
        "--only-mapped-target",
    )

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\tZ\tBETA\tSE\tDirection\tEAF\n"
        "1\t100\trs1\tG\tA\t0.05\t-1.95996398454005\t-0.5\t0.255106728462327\t+-\t0.8\n"
    )

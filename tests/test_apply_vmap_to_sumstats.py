import argparse
import gzip
import math

import pandas as pd
import yaml

from genomatch.apply_vmap_sumstats import run_sumstats_apply
from genomatch.sumstats_utils import iter_sumstats_table_chunks

from utils import REPO_ROOT, read_tsv, run_py, write_json, write_lines


def write_vmap_with_meta(path, lines):
    write_lines(path, lines)
    write_json(
        path.with_name(path.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )


def read_yaml(path):
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def assert_raw_sumstats_metadata_schema_shape(metadata):
    schema = read_yaml(REPO_ROOT / "src/genomatch/schemas/raw-sumstats-metadata.yaml")
    properties = schema["properties"]
    assert set(metadata).issubset(properties)
    for key in schema["required"]:
        assert key in metadata
    for key, value in metadata.items():
        property_schema = properties[key]
        if "enum" in property_schema:
            assert value in property_schema["enum"]
        if property_schema.get("type") == "string":
            assert isinstance(value, str)
        elif property_schema.get("type") == "boolean":
            assert isinstance(value, bool)
    assert "col_SNP" in metadata or {"col_CHR", "col_POS"}.issubset(metadata)
    assert any(key in metadata for key in (
        "col_BETA",
        "col_SE",
        "col_Z",
        "col_P",
        "col_OR",
        "col_ORL95",
        "col_ORU95",
        "col_N",
        "col_CaseN",
        "col_ControlN",
        "col_StudyN",
        "col_EAF",
        "col_OAF",
        "col_CaseEAF",
        "col_CaseOAF",
        "col_ControlEAF",
        "col_ControlOAF",
        "col_INFO",
        "col_Direction",
    ))


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
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\tOR\tEAF\n"
        "1\t200\t1:200:T:C\tT\tC\t-1\t0.25\t0.7\n"
        "1\t100\t1:100:A:G\tA\tG\t0.5\t2.0\t0.1\n"
    )
    sidecar = read_yaml(out.with_name(out.name + ".meta.yaml"))
    assert_raw_sumstats_metadata_schema_shape(sidecar)
    assert sidecar == {
        "cleansumstats_metafile_kind": "minimal",
        "path_sumStats": "out.tsv",
        "delimiter": "\t",
        "missing_value": "",
        "genome_build": "GRCh37",
        "contig_naming": "ncbi",
        "stats_Model": "other",
        "col_CHR": "CHR",
        "col_POS": "POS",
        "col_SNP": "SNP",
        "col_EffectAllele": "EffectAllele",
        "col_OtherAllele": "OtherAllele",
        "clean": False,
        "col_BETA": "BETA",
        "col_OR": "OR",
        "col_EAF": "EAF",
    }


def test_apply_vmap_to_sumstats_projection_drops_unrecognized_input_columns(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA\tUNUSED", "1\t100\tA\tG\t0.5\tdropme"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_vmap_with_meta(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n1\t100\t1:100:A:G\tA\tG\t0.5\n"


def test_apply_vmap_to_sumstats_projection_applies_metadata_missing_value(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA\tG\t-9"])
    write_lines(
        meta,
        [
            'missing_value: "-9"',
            "col_CHR: CHR",
            "col_POS: POS",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_BETA: BETA",
        ],
    )
    write_vmap_with_meta(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n1\t100\t1:100:A:G\tA\tG\t\n"


def test_apply_vmap_to_sumstats_projection_normalizes_common_missing_tokens(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(
        sumstats,
        [
            "CHR\tPOS\tEA\tOA\tBETA",
            "1\t100\tA\tG\tNA",
            "1\t200\tA\tG\tN/A",
            "1\t300\tA\tG\tnull",
            "1\t400\tA\tG\t.",
        ],
    )
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_vmap_with_meta(
        vmap,
        [
            "1\t100\tt1\tA\tG\t.\t0\tidentity",
            "1\t200\tt2\tA\tG\t.\t1\tidentity",
            "1\t300\tt3\tA\tG\t.\t2\tidentity",
            "1\t400\tt4\tA\tG\t.\t3\tidentity",
        ],
    )

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n"
        "1\t100\t1:100:A:G\tA\tG\t\n"
        "1\t200\t1:200:A:G\tA\tG\t\n"
        "1\t300\t1:300:A:G\tA\tG\t\n"
        "1\t400\t1:400:A:G\tA\tG\t\n"
    )


def test_apply_vmap_to_sumstats_projection_rejects_ambiguous_normalized_header_lookup(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA\t beta ", "1\t100\tA\tG\t0.5\t0.6"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: beta"])
    write_vmap_with_meta(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)

    assert result.returncode != 0
    assert "identified uniquely" in result.stderr


def test_apply_vmap_to_sumstats_projection_rejects_duplicate_payload_output_names(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA\tG\t0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA", "col_Z: BETA"])
    write_vmap_with_meta(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)

    assert result.returncode != 0
    assert "duplicate output name" in result.stderr


def test_apply_vmap_to_sumstats_projection_rejects_payload_collision_with_canonical_variant_column(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA\tG\t0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: POS"])
    write_vmap_with_meta(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity"])

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)

    assert result.returncode != 0
    assert "would collide with canonical output column" in result.stderr


def test_apply_vmap_to_sumstats_resolves_padded_headers_for_projection_output(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR \t POS \t EA \t OA \t BETA ", "1\t100\tA\tG\t0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(vmap, ["2\t250\tt1\tT\tC\t.\t0\tswap"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n2\t250\t2:250:T:C\tT\tC\t-0.5\n"


def test_apply_vmap_to_sumstats_uses_canonical_variant_columns_even_when_source_variant_columns_are_missing(tmp_path):
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
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n1\t100\t1:100:A:G\tA\tG\t0.5\n"


def test_apply_vmap_to_sumstats_rejects_missing_provenance_sentinel_rows(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA\tG\t0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(vmap, ["1\t100\tt1\tA\tG\t.\t0\tidentity", "1\t200\tt2\tC\tT\t.\t-1\tmissing"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode != 0
    assert "vmap row source_index out of range" in result.stderr


def test_apply_vmap_to_sumstats_warns_and_writes_empty_for_invalid_swapped_numeric_effects(tmp_path):
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
    assert "writing missing value" in result.stderr.lower()
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\tOR\tEAF\tORL95\tORU95\n"
        "1\t100\t1:100:G:A\tG\tA\t\t\t\t\t\n"
    )


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
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\tOR\tEAF\n1\t100\t1:100:G:T\tG\tT\t-0.5\t0.5\t0.9\n"


def test_apply_vmap_to_sumstats_swap_complements_other_allele_frequency_columns(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tOAF\tCaseOAF\tControlOAF", "1\t100\tA\tG\t0.2\t0.3\t0.4"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_OAF: OAF",
            "col_CaseOAF: CaseOAF",
            "col_ControlOAF: ControlOAF",
        ],
    )
    write_lines(vmap, ["1\t100\tt1\tG\tA\t.\t0\tswap"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tOAF\tCaseOAF\tControlOAF\n"
        "1\t100\t1:100:G:A\tG\tA\t0.8\t0.7\t0.6\n"
    )


def test_apply_vmap_to_sumstats_formats_swapped_object_floats_with_float_format(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(
        sumstats,
        [
            "CHR\tPOS\tEA\tOA\tBETA",
            "1\t100\tA\tG\t0.12345678901234567",
            "1\t200\tC\tT\t0.12345678901234567",
        ],
    )
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(
        vmap,
        [
            "1\t100\tt1\tA\tG\t.\t0\tidentity",
            "1\t200\tt2\tT\tC\t.\t1\tswap",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n"
        "1\t100\t1:100:A:G\tA\tG\t0.12345678901234567\n"
        "1\t200\t1:200:T:C\tT\tC\t-0.123456789012346\n"
    )


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
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n2\t250\t2:250:T:C\tT\tC\t-0.5\n"


def test_apply_vmap_to_sumstats_retain_snp_id_uses_vmap_id(tmp_path):
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
    write_lines(vmap, ["2\t250\trs_target\tT\tC\t.\t0\tidentity"])
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
        "--retain-snp-id",
    )

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n2\t250\trs_target\tT\tC\t0.5\n"


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
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tORL95\tORU95\n"
        "1\t100\t1:100:G:A\tG\tA\t0.666666666666667\t0.833333333333333\n"
    )


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
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n1\t100\t1:100:A:G\tA\tG\t0.5\n"


def test_apply_vmap_to_sumstats_does_not_rewrite_joined_chr_field_in_place(tmp_path):
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
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n2\t100\t2:100:A:G\tA\tG\t0.5\n"


def test_apply_vmap_to_sumstats_writes_tab_delimited_projection_output(tmp_path):
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
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n1\t100\t1:100:G:A\tG\tA\t-0.5\n"


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
        assert handle.read() == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n1\t100\t1:100:G:A\tG\tA\t-0.5\n"


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


def test_apply_vmap_to_sumstats_surfaces_pandas_parser_failure_for_ragged_payload_row(tmp_path):
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
    assert "Too many columns specified" in result.stderr


def test_apply_vmap_to_sumstats_roundtrip_preserves_provenance_with_comment_and_blank_lines(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "imported.vmap"
    out = tmp_path / "out.tsv"

    sumstats.write_text(
        "# metadata preamble\n"
        "\n"
        "CHR\tPOS\tSNP\tEA\tOA\tBETA\n"
        "# in-body comment that should not affect source_index\n"
        "1\t100\trs1\tA\tG\t0.5\n"
        "\n"
        "# another in-body comment\n"
        "1\t200\trs2\tC\tT\t-1.0\n",
        encoding="utf-8",
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
        ],
    )

    result = run_py("import_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--output", vmap)
    assert result.returncode == 0, result.stderr
    assert read_tsv(vmap) == [
        ["1", "100", "rs1", "A", "G", ".", "0", "identity"],
        ["1", "200", "rs2", "C", "T", ".", "1", "identity"],
    ]

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n"
        "1\t100\t1:100:A:G\tA\tG\t0.5\n"
        "1\t200\t1:200:C:T\tC\tT\t-1.0\n"
    )


def test_apply_vmap_to_sumstats_source_index_ignores_comment_and_blank_lines(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"

    sumstats.write_text(
        "# pre-header comment\n"
        "\n"
        "CHR\tPOS\tSNP\tEA\tOA\tBETA\n"
        "# skipped line before first data row\n"
        "1\t100\trs1\tA\tG\t0.5\n"
        "\n"
        "# skipped line before second data row\n"
        "1\t200\trs2\tC\tT\t-1.0\n",
        encoding="utf-8",
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
        ],
    )
    write_lines(vmap, ["1\t200\trs2\tC\tT\t.\t1\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n1\t200\t1:200:C:T\tC\tT\t-1.0\n"


def test_sumstats_chunk_iterator_preserves_source_index_across_chunks(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    sumstats.write_text(
        "# pre-header comment\n"
        "CHR\tPOS\tBETA\n"
        "# skipped comment\n"
        "1\t100\t0.5\n"
        "\n"
        "1\t200\t-1.0\n"
        "1\t300\t2.0\n",
        encoding="utf-8",
    )

    chunks = list(iter_sumstats_table_chunks(sumstats, chunk_rows=1))

    assert [chunk.source_index.tolist() for chunk in chunks] == [[0], [1], [2]]
    assert [chunk.frame["POS"].tolist() for chunk in chunks] == [["100"], ["200"], ["300"]]


def test_run_sumstats_apply_finds_vmap_rows_after_first_chunk(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    out = tmp_path / "out.tsv"
    write_lines(
        sumstats,
        [
            "CHR\tPOS\tEA\tOA\tBETA",
            "1\t100\tA\tG\t0.1",
            "1\t200\tA\tG\t0.2",
            "1\t300\tA\tG\t0.3",
            "1\t400\tA\tG\t0.4",
        ],
    )
    metadata = {
        "col_CHR": "CHR",
        "col_POS": "POS",
        "col_EffectAllele": "EA",
        "col_OtherAllele": "OA",
        "col_BETA": "BETA",
    }
    vmap_frame = pd.DataFrame(
        {
            "chrom": ["1", "1"],
            "pos": ["300", "400"],
            "id": ["t3", "t4"],
            "a1": ["A", "G"],
            "a2": ["G", "A"],
            "source_shard": [".", "."],
            "source_index": [2, 3],
            "allele_op": ["identity", "swap"],
        }
    )

    result = run_sumstats_apply(
        argparse.Namespace(clean=False, fill_mode="column", use_af_inference=False, retain_snp_id=False),
        metadata=metadata,
        retained_columns=[("col_BETA", 4, "BETA")],
        projection_payload_mappings={"col_BETA": "BETA"},
        vmap_frame=vmap_frame,
        vmap_meta={"target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
        input_path=sumstats,
        input_header=["CHR", "POS", "EA", "OA", "BETA"],
        output_path=out,
        chunk_rows=2,
    )

    assert result == 0
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n"
        "1\t300\t1:300:A:G\tA\tG\t0.3\n"
        "1\t400\t1:400:G:A\tG\tA\t-0.4\n"
    )


def test_apply_vmap_to_sumstats_reuses_one_source_row_for_multiple_vmap_rows(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA\tBETA", "1\t100\tA\tG\t0.5"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA", "col_BETA: BETA"])
    write_lines(
        vmap,
        [
            "1\t100\tt1\tA\tG\t.\t0\tidentity",
            "2\t250\tt2\tT\tC\t.\t0\tswap",
        ],
    )
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tBETA\n"
        "1\t100\t1:100:A:G\tA\tG\t0.5\n"
        "2\t250\t2:250:T:C\tT\tC\t-0.5\n"
    )


def test_apply_vmap_to_sumstats_supports_hash_prefixed_chr_header(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(
        sumstats,
        [
            "#chrom\tpos\tref\talt\trsids\tbeta",
            "1\t13668\tG\tA\trs2691328\t0.2",
        ],
    )
    write_lines(
        meta,
        [
            'col_CHR: "#chrom"',
            "col_POS: pos",
            "col_SNP: rsids",
            "col_EffectAllele: alt",
            "col_OtherAllele: ref",
            "col_BETA: beta",
        ],
    )
    write_lines(vmap, ["1\t13668\trs2691328\tA\tG\t.\t0\tidentity"])
    write_json(vmap.with_name(vmap.name + ".meta.json"), {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}})

    result = run_py("apply_vmap_to_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--vmap", vmap, "--output", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tbeta\n1\t13668\t1:13668:A:G\tA\tG\t0.2\n"


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
        "1\t100\t1:100:T:C\tT\tC\t0.05\t1.95996398454005\t0.5\t0.255106728462327\n"
    )
    sidecar = read_yaml(out.with_name(out.name + ".meta.yaml"))
    assert_raw_sumstats_metadata_schema_shape(sidecar)
    assert sidecar["path_sumStats"] == "out.tsv"
    assert sidecar["delimiter"] == "\t"
    assert sidecar["missing_value"] == ""
    assert sidecar["genome_build"] == "GRCh37"
    assert sidecar["contig_naming"] == "ncbi"
    assert sidecar["clean"] is True
    assert sidecar["fill_mode"] == "column"
    assert sidecar["use_af_inference"] is False
    assert sidecar["col_P"] == "P"
    assert sidecar["col_BETA"] == "BETA"
    assert "payload_mappings" not in sidecar


def test_apply_vmap_to_sumstats_clean_rejects_duplicate_payload_output_names_before_harmonization(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t100\trs1\tA\tG\t0.5"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_BETA: BETA",
            "col_Z: BETA",
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
    assert "duplicate output name" in result.stderr


def test_apply_vmap_to_sumstats_clean_rejects_payload_collision_with_canonical_variant_column(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(sumstats, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t100\trs1\tA\tG\t0.5"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_BETA: SNP",
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
    assert "would collide with canonical output column" in result.stderr


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
    assert out_log.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\n1\t100\t1:100:A:G\tA\tG\t\n"


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
        "1\t100\t1:100:A:G\tA\tG\t\t\t2\t1",
        "1\t200\t1:200:C:T\tC\tT\t2.25717681190766e-19\t9\t4\t2",
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
        "1\t100\t1:100:A:G\tA\tG\t0.0455002638963584\t2\t2\t1",
        "1\t200\t1:200:C:T\tC\tT\t2.25717681190766e-19\t9\t4\t2",
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
        "1\t100\t1:100:A:G\tA\tG\t0.8\n"
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


def test_apply_vmap_to_sumstats_clean_logistic_vectorizes_or_and_interval_derivations(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(
        sumstats,
        [
            "CHR\tPOS\tSNP\tEA\tOA\tOR\tORL95\tORU95",
            "1\t100\trs1\tA\tG\t2\t1.2\t1.5",
            "1\t200\trs2\tC\tT\t4\t3\t5",
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
            "col_OR: OR",
            "col_ORL95: ORL95",
            "col_ORU95: ORU95",
            "stats_Model: logistic",
        ],
    )
    write_vmap_with_meta(
        vmap,
        [
            "1\t100\trs1\tA\tG\t.\t0\tidentity",
            "1\t200\trs2\tC\tT\t.\t1\tidentity",
        ],
    )

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
    lines = out.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    row1 = dict(zip(header, lines[1].split("\t")))
    row2 = dict(zip(header, lines[2].split("\t")))
    assert math.isclose(float(row1["BETA"]), math.log(2.0), rel_tol=1e-12)
    assert math.isclose(float(row2["BETA"]), math.log(4.0), rel_tol=1e-12)
    assert math.isclose(float(row1["SE"]), 0.056925421353233446, rel_tol=1e-12)
    assert math.isclose(float(row2["SE"]), 0.13031505369366933, rel_tol=1e-12)


def test_apply_vmap_to_sumstats_clean_use_af_inference_vectorizes_beta_and_se(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    vmap = tmp_path / "map.vmap"
    out = tmp_path / "out.tsv"
    write_lines(
        sumstats,
        [
            "CHR\tPOS\tSNP\tEA\tOA\tZ\tN\tEAF",
            "1\t100\trs1\tA\tG\t2\t100\t0.25",
            "1\t200\trs2\tC\tT\t3\t200\t0.4",
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
            "col_Z: Z",
            "col_N: N",
            "col_EAF: EAF",
            "stats_Model: linear",
        ],
    )
    write_vmap_with_meta(
        vmap,
        [
            "1\t100\trs1\tA\tG\t.\t0\tidentity",
            "1\t200\trs2\tC\tT\t.\t1\tidentity",
        ],
    )

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
    lines = out.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    row1 = dict(zip(header, lines[1].split("\t")))
    row2 = dict(zip(header, lines[2].split("\t")))
    assert math.isclose(float(row1["BETA"]), 0.32025630761017426, rel_tol=1e-12)
    assert math.isclose(float(row2["BETA"]), 0.29952114893657694, rel_tol=1e-12)
    assert math.isclose(float(row1["SE"]), 0.16012815380508713, rel_tol=1e-12)
    assert math.isclose(float(row2["SE"]), 0.09984038297885897, rel_tol=1e-12)


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


def test_apply_vmap_to_sumstats_clean_fills_metadata_total_n_for_mapped_row(tmp_path):
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
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tN\tBETA\n"
        "1\t100\t1:100:A:G\tA\tG\t100\t0.5\n"
    )


def test_apply_vmap_to_sumstats_clean_swap_keeps_direction(tmp_path):
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
    )

    assert result.returncode == 0, result.stderr
    assert out.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\tZ\tBETA\tSE\tDirection\tEAF\n"
        "1\t100\t1:100:G:A\tG\tA\t0.05\t-1.95996398454005\t-0.5\t0.255106728462327\t+-\t0.8\n"
        "1\t200\t1:200:T:C\tT\tC\t\t\t1.5\t\t--\t0.4\n"
    )

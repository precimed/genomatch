import json
import pytest

from utils import run_py_with_env, write_fasta, write_json, write_lines, write_match_config


def test_guess_build_updates_unknown_build_only_for_vtable(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "AC"})
    write_fasta(grch38, {"chr1": "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\trs1\tA\tC", "1\t2\trs2\tC\tA"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ncbi"},
    )

    result = run_py_with_env("guess_build.py", {"MATCH_CONFIG": str(config)}, "--input", source, "--write")
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary["genome_build"] == "GRCh37"
    meta = json.loads(source.with_name(source.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta["genome_build"] == "GRCh37"


def test_guess_build_updates_vmap_target_metadata_only(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vmap"
    write_fasta(grch37, {"chr1": "AC"})
    write_fasta(grch38, {"chr1": "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\trs1\tA\tC\tshardA\t0\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {
            "object_type": "variant_map",
            "target": {"genome_build": "unknown", "contig_naming": "ncbi"},
            "extra": "kept",
        },
    )

    result = run_py_with_env("guess_build.py", {"MATCH_CONFIG": str(config)}, "--input", source, "--write")
    assert result.returncode == 0, result.stderr
    meta = json.loads(source.with_name(source.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta["target"]["genome_build"] == "GRCh37"
    assert meta["target"]["contig_naming"] == "ncbi"
    assert meta["extra"] == "kept"
    read_text = source.read_text(encoding="utf-8")
    assert "shardA\t0\tidentity" in read_text


def test_guess_build_force_overwrites_existing_build(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "AC"})
    write_fasta(grch38, {"chr1": "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\trs1\tA\tC", "1\t2\trs2\tC\tA"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": "ncbi"},
    )

    result = run_py_with_env("guess_build.py", {"MATCH_CONFIG": str(config)}, "--input", source, "--write", "--force")
    assert result.returncode == 0, result.stderr
    meta = json.loads(source.with_name(source.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta["genome_build"] == "GRCh37"


def test_guess_build_requires_declared_contig_naming(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "AC"})
    write_fasta(grch38, {"chr1": "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\trs1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown"},
    )

    result = run_py_with_env("guess_build.py", {"MATCH_CONFIG": str(config)}, "--input", source)
    assert result.returncode != 0
    assert "normalize_contigs.py" in result.stderr


def test_guess_build_returns_unknown_on_ambiguous_build_evidence(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "A"})
    write_fasta(grch38, {"chr1": "A"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\trs1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ncbi"},
    )

    result = run_py_with_env("guess_build.py", {"MATCH_CONFIG": str(config)}, "--input", source)
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary["genome_build"] == "unknown"
    assert summary["confidence"] == "low"


def test_guess_build_accepts_ucsc_named_input(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "AC"})
    write_fasta(grch38, {"chr1": "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["chr1\t1\trs1\tA\tC", "chr1\t2\trs2\tC\tA"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ucsc"},
    )

    result = run_py_with_env("guess_build.py", {"MATCH_CONFIG": str(config)}, "--input", source)
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary["genome_build"] == "GRCh37"
    assert summary["contig_naming"] == "ucsc"
    assert summary["normalization"] == "none"


def test_guess_build_returns_unknown_when_best_rate_is_below_threshold(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "AC"})
    write_fasta(grch38, {"chr1": "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\trs1\tA\tC", "1\t2\trs2\tA\tA"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ncbi"},
    )

    result = run_py_with_env("guess_build.py", {"MATCH_CONFIG": str(config)}, "--input", source)
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary["genome_build"] == "unknown"
    assert summary["confidence"] == "low"


def test_guess_build_returns_unknown_when_no_rows_are_usable_as_evidence(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "AC"})
    write_fasta(grch38, {"chr1": "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t999\trs1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ncbi"},
    )

    result = run_py_with_env("guess_build.py", {"MATCH_CONFIG": str(config)}, "--input", source)
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary["genome_build"] == "unknown"
    assert summary["confidence"] == "low"
    assert [entry["considered_rows"] for entry in summary["evidence"]] == [0, 0]


def test_guess_build_fails_on_malformed_sidecar(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "AC"})
    write_fasta(grch38, {"chr1": "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\trs1\tA\tC"])
    source.with_name(source.name + ".meta.json").write_text("{bad json\n", encoding="utf-8")

    result = run_py_with_env("guess_build.py", {"MATCH_CONFIG": str(config)}, "--input", source, "--write")
    assert result.returncode != 0
    assert "json" in result.stderr.lower() or "expecting property name enclosed in double quotes" in result.stderr.lower()


@pytest.mark.parametrize(
    "contig_naming,chrom,pos",
    [
        ("plink", "23", "60001"),
        ("plink_splitx", "25", "60001"),
    ],
)
def test_guess_build_accepts_plink_family_input(tmp_path, contig_naming, chrom, pos):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / f"cohort.{contig_naming}.vtable"
    prefix = "N" * 60000
    write_fasta(grch37, {"chrX": prefix + "AC"})
    write_fasta(grch38, {"chrX": prefix + "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, [f"{chrom}\t{pos}\trs1\tA\tC", f"{chrom}\t60002\trs2\tC\tA"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": contig_naming},
    )

    result = run_py_with_env("guess_build.py", {"MATCH_CONFIG": str(config)}, "--input", source)
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary["genome_build"] == "GRCh37"
    assert summary["contig_naming"] == contig_naming
    assert summary["normalization"] == f"{contig_naming}->ucsc"


def test_guess_build_reference_access_modes_are_consistent_and_default_to_bulk(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "ACG"})
    write_fasta(grch38, {"chr1": "TGA"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\trs1\tA\tC", "1\t2\trs2\tC\tA", "1\t3\trs3\tT\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ncbi"},
    )

    common_env = {"MATCH_CONFIG": str(config)}
    result_default = run_py_with_env("guess_build.py", common_env, "--input", source)
    result_bulk = run_py_with_env(
        "guess_build.py",
        {**common_env, "MATCH_REFERENCE_ACCESS_MODE": "BULK"},
        "--input",
        source,
    )
    result_legacy = run_py_with_env(
        "guess_build.py",
        {**common_env, "MATCH_REFERENCE_ACCESS_MODE": "legacy"},
        "--input",
        source,
    )
    assert result_default.returncode == 0, result_default.stderr
    assert result_bulk.returncode == 0, result_bulk.stderr
    assert result_legacy.returncode == 0, result_legacy.stderr

    summary_default = json.loads(result_default.stdout)
    summary_bulk = json.loads(result_bulk.stdout)
    summary_legacy = json.loads(result_legacy.stdout)
    assert summary_default == summary_bulk
    assert summary_bulk == summary_legacy


def test_guess_build_sample_rows_limits_considered_rows(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "A" * 32})
    write_fasta(grch38, {"chr1": "T" * 32})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    rows = [f"1\t{i}\trs{i}\tA\tC" for i in range(1, 21)]
    write_lines(source, rows)
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "guess_build.py",
        {"MATCH_CONFIG": str(config)},
        "--input",
        source,
        "--sample-rows",
        "5",
    )
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary["sample_rows_requested"] == 5
    assert summary["sample_rows_used"] == 5
    considered = [entry["considered_rows"] for entry in summary["evidence"]]
    assert considered == [5, 5]


def test_guess_build_sample_rows_zero_uses_all_rows(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "A" * 32})
    write_fasta(grch38, {"chr1": "T" * 32})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    rows = [f"1\t{i}\trs{i}\tA\tC" for i in range(1, 13)]
    write_lines(source, rows)
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "guess_build.py",
        {"MATCH_CONFIG": str(config)},
        "--input",
        source,
        "--sample-rows",
        "0",
    )
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary["sample_rows_requested"] == 0
    assert summary["sample_rows_used"] == 12
    considered = [entry["considered_rows"] for entry in summary["evidence"]]
    assert considered == [12, 12]


def test_guess_build_rejects_negative_sample_rows(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "cohort.vtable"
    write_fasta(grch37, {"chr1": "AC"})
    write_fasta(grch38, {"chr1": "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\trs1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ncbi"},
    )

    result = run_py_with_env("guess_build.py", {"MATCH_CONFIG": str(config)}, "--input", source, "--sample-rows", "-1")
    assert result.returncode != 0
    assert "--sample-rows must be >= 0" in result.stderr

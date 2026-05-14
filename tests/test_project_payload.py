from pathlib import Path

import numpy as np
import pgenlib

from utils import read_bed_matrix, read_bim, run_py, write_bfile, write_json, write_lines


def assert_wrote(result, value: Path) -> None:
    assert f"project_payload.py: wrote {value}" in result.stderr


def write_sumstats_metadata(path: Path) -> None:
    write_lines(
        path,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_BETA: BETA",
        ],
    )


def write_vmap(path: Path, lines: list[str], *, genome_build: str = "GRCh38", contig_naming: str = "ncbi") -> None:
    write_lines(path, lines)
    write_json(
        path.with_name(path.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": genome_build, "contig_naming": contig_naming}},
    )


def write_psam(path: Path, sample_lines: list[str], *, header: str = "#IID") -> None:
    write_lines(path, [header, *sample_lines])


def write_pvar(path: Path, variant_lines: list[str]) -> None:
    write_lines(path, ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT", *variant_lines])


def write_pfile(
    prefix: Path,
    variant_lines: list[str],
    sample_lines: list[str],
    rows: list[dict],
    *,
    psam_header: str = "#IID",
) -> None:
    write_pvar(prefix.with_suffix(".pvar"), variant_lines)
    write_psam(prefix.with_suffix(".psam"), sample_lines, header=psam_header)
    writer = pgenlib.PgenWriter(
        str(prefix.with_suffix(".pgen")).encode("utf-8"),
        len(sample_lines),
        variant_ct=len(rows),
        nonref_flags=False,
        allele_ct_limit=2,
        hardcall_phase_present=any(row["channel"] == "phase" for row in rows),
        dosage_present=any(row["channel"] == "dosage" for row in rows),
    )
    try:
        for row in rows:
            if row["channel"] == "dosage":
                writer.append_dosages(np.array(row["dosages"], dtype=np.float32))
            elif row["channel"] == "phase":
                writer.append_partially_phased(
                    np.array(row["alleles"], dtype=np.int32),
                    np.array(row["phasepresent"], dtype=np.uint8),
                )
            else:
                writer.append_alleles(np.array(row["alleles"], dtype=np.int32))
    finally:
        writer.close()


def read_pfile_genotypes(prefix: Path, n_samples: int, variant_idx: int) -> list[int]:
    reader = pgenlib.PgenReader(str(prefix.with_suffix(".pgen")).encode("utf-8"), n_samples)
    try:
        out = np.empty(n_samples, dtype=np.int8)
        reader.read(variant_idx, out)
        return out.astype(int).tolist()
    finally:
        reader.close()


def test_project_payload_requires_vmap(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    output = tmp_path / "aligned.tsv"
    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs1\tG\tA\t0.5"])
    write_sumstats_metadata(metadata)

    result = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--output",
        output,
    )

    assert result.returncode != 0
    assert "the following arguments are required: --vmap" in result.stderr


def test_project_payload_rejects_legacy_mapping_flags(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"
    output = tmp_path / "aligned.tsv"
    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs1\tG\tA\t0.5"])
    write_sumstats_metadata(metadata)
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])

    base = [
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        output,
    ]
    cases = [
        ("--target", tmp_path / "target.vtable"),
        ("--source-vmap", vmap),
        ("--full-target", None),
    ]

    for flag, value in cases:
        args = [*base]
        args.append(flag)
        if value is not None:
            args.append(value)
        result = run_py(*args)
        assert result.returncode != 0
        assert f"unrecognized arguments: {flag}" in result.stderr


def test_project_payload_sumstats_uses_explicit_vmap_and_does_not_retain_restricted_vmap(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"
    output = tmp_path / "aligned.tsv"

    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs1\tG\tA\t0.5"])
    write_sumstats_metadata(metadata)
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])

    result = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert "match_vmap_to_target.py" not in result.stderr
    assert "restrict_vmap.py" not in result.stderr
    assert "apply_vmap_to_sumstats.py" in result.stderr
    assert str(vmap) in result.stderr
    assert not (tmp_path / "aligned.tsv.vmap").exists()
    assert output.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEA\tOA\tBETA\n1\t1\t1:1:G:A\tG\tA\t0.5\n"
    assert_wrote(result, output)


def test_project_payload_retain_snp_id_passes_through_to_apply(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"
    output = tmp_path / "aligned.tsv"

    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs_source\tG\tA\t0.5"])
    write_sumstats_metadata(metadata)
    write_vmap(vmap, ["1\t1\tretained_target_id\tG\tA\t.\t0\tidentity"])

    result = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        output,
        "--retain-snp-id",
    )

    assert result.returncode == 0, result.stderr
    assert "--retain-snp-id" in result.stderr
    assert output.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEA\tOA\tBETA\n1\t1\tretained_target_id\tG\tA\t0.5\n"


def test_project_payload_sumstats_clean_dispatches_clean_apply(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"
    output = tmp_path / "aligned.tsv"

    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tP\tBETA", "1\t1\trs1\tG\tA\t0.05\t0.5"])
    write_lines(
        metadata,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_P: P",
            "col_BETA: BETA",
            "stats_Model: linear",
        ],
    )
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])

    result = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats-clean",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert "apply_vmap_to_sumstats.py" in result.stderr
    assert "--clean" in result.stderr
    assert output.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\tZ\tBETA\tSE\n"
        "1\t1\t1:1:G:A\tG\tA\t0.05\t1.95996398454005\t0.5\t0.255106728462327\n"
    )


def test_project_payload_sumstats_uses_metadata_path_sumstats_when_input_omitted(tmp_path):
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    source = bundle / "study.tsv"
    metadata = bundle / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"
    output = tmp_path / "aligned.tsv"

    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs1\tG\tA\t0.5"])
    write_lines(
        metadata,
        [
            "path_sumStats: study.tsv",
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_BETA: BETA",
        ],
    )
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])

    result = run_py(
        "project_payload.py",
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    apply_line = next(line for line in result.stderr.splitlines() if "apply_vmap_to_sumstats.py" in line)
    assert "--input" not in apply_line
    assert output.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEA\tOA\tBETA\n1\t1\t1:1:G:A\tG\tA\t0.5\n"


def test_project_payload_sumstats_clean_passes_fill_mode_and_af_inference(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"
    output = tmp_path / "aligned.tsv"

    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tZ\tBETA\tEAF", "1\t1\trs1\tG\tA\t2.0\t0.5\t0.4"])
    write_lines(
        metadata,
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
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])

    result = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats-clean",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        output,
        "--fill-mode",
        "row",
        "--use-af-inference",
    )

    assert result.returncode == 0, result.stderr
    assert "--clean" in result.stderr
    assert "--fill-mode row" in result.stderr
    assert "--use-af-inference" in result.stderr
    assert output.read_text(encoding="utf-8") == (
        "CHR\tPOS\tSNP\tEffectAllele\tOtherAllele\tP\tZ\tN\tBETA\tSE\tEAF\n"
        "1\t1\t1:1:G:A\tG\tA\t0.0455002638963584\t2\t29.3333333333333\t0.5\t0.25\t0.4\n"
    )


def test_project_payload_bfile_uses_explicit_vmap_and_allows_sharded_output_prefix(tmp_path):
    source_prefix = tmp_path / "source"
    source_bim = source_prefix.with_suffix(".bim")
    vmap = tmp_path / "cohort.prepared.vmap"
    output = tmp_path / "aligned.@"

    write_bfile(
        source_prefix,
        ["1\trs1\t0\t1\tG\tA", "2\trs2\t0\t2\tC\tA"],
        ["S1 S1 0 0 0 -9"],
        [[0], [3]],
    )
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity", "2\t2\trs2\tC\tA\t.\t1\tidentity"])

    result = run_py(
        "project_payload.py",
        "--input",
        source_bim,
        "--input-format",
        "bfile",
        "--vmap",
        vmap,
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert "apply_vmap_to_bfile.py" in result.stderr
    assert "match_vmap_to_target.py" not in result.stderr
    assert not (tmp_path / "aligned.all_targets.vmap").exists()
    assert (tmp_path / "aligned.1.bim").exists()
    assert (tmp_path / "aligned.2.bim").exists()
    assert [row.snp for row in read_bim(tmp_path / "aligned.1.bim")] == ["1:1:G:A"]
    assert [row.snp for row in read_bim(tmp_path / "aligned.2.bim")] == ["2:2:C:A"]
    assert read_bed_matrix(tmp_path / "aligned.1.bed", 1, 1) == [[0]]
    assert read_bed_matrix(tmp_path / "aligned.2.bed", 1, 1) == [[3]]
    assert_wrote(result, output)


def test_project_payload_pfile_runs_apply_step_and_writes_projected_payload(tmp_path):
    source_prefix = tmp_path / "source"
    source_pvar = source_prefix.with_suffix(".pvar")
    vmap = tmp_path / "cohort.prepared.vmap"
    output = tmp_path / "aligned"

    write_pfile(
        source_prefix,
        ["1\t1\trs1\tA\tG"],
        ["S1", "S2"],
        [{"channel": "hardcall", "alleles": [0, 0, 0, 1]}],
    )
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])

    result = run_py(
        "project_payload.py",
        "--input",
        source_pvar,
        "--input-format",
        "pfile",
        "--vmap",
        vmap,
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert "apply_vmap_to_pfile.py" in result.stderr
    assert not (tmp_path / "aligned.vmap").exists()
    assert output.with_suffix(".pgen").exists()
    assert output.with_suffix(".pvar").exists()
    assert output.with_suffix(".psam").exists()
    assert read_pfile_genotypes(output, 2, 0) == [0, 1]
    assert_wrote(result, output)


def test_project_payload_validates_common_argument_combinations(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"

    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs1\tG\tA\t0.5"])
    write_sumstats_metadata(metadata)
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])

    missing_metadata = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--vmap",
        vmap,
        "--output",
        tmp_path / "aligned.tsv",
    )
    assert missing_metadata.returncode != 0
    assert "--sumstats-metadata is required for --input-format=sumstats" in missing_metadata.stderr

    rejected_output = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        tmp_path / "aligned.@.tsv",
    )
    assert rejected_output.returncode != 0
    assert "--output must not contain @ for --input-format=sumstats" in rejected_output.stderr

    fill_mode_rejected = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        tmp_path / "aligned-fill.tsv",
        "--fill-mode",
        "row",
    )
    assert fill_mode_rejected.returncode != 0
    assert "--fill-mode is supported only for --input-format=sumstats-clean" in fill_mode_rejected.stderr

    af_rejected = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        tmp_path / "aligned-af.tsv",
        "--use-af-inference",
    )
    assert af_rejected.returncode != 0
    assert "--use-af-inference is supported only for --input-format=sumstats-clean" in af_rejected.stderr


def test_project_payload_rejects_non_mapped_only_vmap(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"
    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs1\tG\tA\t0.5"])
    write_sumstats_metadata(metadata)
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t-1\tmissing"])

    result = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        tmp_path / "aligned.tsv",
    )

    assert result.returncode != 0
    assert "is not a mapped-only .vmap" in result.stderr


def test_project_payload_rejects_genotype_metadata_and_wrong_target_sample_flags(tmp_path):
    source_prefix = tmp_path / "source"
    source_bim = source_prefix.with_suffix(".bim")
    vmap = tmp_path / "cohort.prepared.vmap"
    metadata = tmp_path / "study.yaml"
    target_psam = tmp_path / "target.psam"

    write_bfile(source_prefix, ["1\trs1\t0\t1\tG\tA"], ["S1 S1 0 0 0 -9"], [[0]])
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])
    write_sumstats_metadata(metadata)
    write_psam(target_psam, ["S1"])

    metadata_rejected = run_py(
        "project_payload.py",
        "--input",
        source_bim,
        "--input-format",
        "bfile",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        tmp_path / "aligned",
    )
    assert metadata_rejected.returncode != 0
    assert "--sumstats-metadata is supported only for --input-format=sumstats or sumstats-clean" in metadata_rejected.stderr

    target_psam_rejected = run_py(
        "project_payload.py",
        "--input",
        source_bim,
        "--input-format",
        "bfile",
        "--target-psam",
        target_psam,
        "--vmap",
        vmap,
        "--output",
        tmp_path / "aligned",
    )
    assert target_psam_rejected.returncode != 0
    assert "--target-psam is supported only for --input-format=pfile" in target_psam_rejected.stderr


def test_project_payload_rejects_resume(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"
    output = tmp_path / "aligned.tsv"

    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs1\tG\tA\t0.5"])
    write_sumstats_metadata(metadata)
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])

    result = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        output,
        "--resume",
    )
    assert result.returncode != 0
    assert "unrecognized arguments: --resume" in result.stderr


def test_project_payload_preflight_excludes_retained_vmap_outputs(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"
    output = tmp_path / "aligned.tsv"

    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs1\tG\tA\t0.5"])
    write_sumstats_metadata(metadata)
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])
    write_vmap(tmp_path / "aligned.tsv.vmap", ["1\t1\trs1\tG\tA\t.\t0\tidentity"])

    result = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "aligned.tsv.vmap").exists()
    assert output.exists()


def test_project_payload_fails_when_wrapper_managed_outputs_exist_without_force(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"
    output = tmp_path / "aligned.tsv"

    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs1\tG\tA\t0.5"])
    write_sumstats_metadata(metadata)
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])
    write_lines(output, ["stale"])

    result = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        output,
    )

    assert result.returncode != 0
    assert "planned output files already exist" in result.stderr
    assert str(output) in result.stderr


def test_project_payload_force_deletes_wrapper_managed_outputs_first_and_reruns(tmp_path):
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    vmap = tmp_path / "study.prepared.vmap"
    output = tmp_path / "aligned.tsv"

    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs1\tG\tA\t0.5"])
    write_sumstats_metadata(metadata)
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])
    write_lines(output, ["stale"])

    result = run_py(
        "project_payload.py",
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--vmap",
        vmap,
        "--output",
        output,
        "--force",
    )

    assert result.returncode == 0, result.stderr
    assert output.read_text(encoding="utf-8") == "CHR\tPOS\tSNP\tEA\tOA\tBETA\n1\t1\t1:1:G:A\tG\tA\t0.5\n"


def test_project_payload_sharded_bfile_force_deletes_only_target_contig_outputs(tmp_path):
    source_prefix = tmp_path / "source"
    source_bim = source_prefix.with_suffix(".bim")
    vmap = tmp_path / "cohort.prepared.vmap"
    output = tmp_path / "aligned.@"

    write_bfile(
        source_prefix,
        ["1\trs1\t0\t1\tG\tA", "2\trs2\t0\t2\tC\tA"],
        ["S1 S1 0 0 0 -9"],
        [[0], [3]],
    )
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity", "2\t2\trs2\tC\tA\t.\t1\tidentity"])
    write_lines(tmp_path / "aligned.1.bim", ["stale"])
    write_lines(tmp_path / "aligned.1.bed", ["stale"])
    write_lines(tmp_path / "aligned.1.fam", ["stale"])
    write_lines(tmp_path / "aligned.2.bim", ["stale"])
    write_lines(tmp_path / "aligned.2.bed", ["stale"])
    write_lines(tmp_path / "aligned.2.fam", ["stale"])
    write_lines(tmp_path / "aligned.unrelated.bim", ["keep"])
    write_lines(tmp_path / "aligned.unrelated.bed", ["keep"])
    write_lines(tmp_path / "aligned.unrelated.fam", ["keep"])

    result = run_py(
        "project_payload.py",
        "--input",
        source_bim,
        "--input-format",
        "bfile",
        "--vmap",
        vmap,
        "--output",
        output,
        "--force",
    )

    assert result.returncode == 0, result.stderr
    assert [row.snp for row in read_bim(tmp_path / "aligned.1.bim")] == ["1:1:G:A"]
    assert [row.snp for row in read_bim(tmp_path / "aligned.2.bim")] == ["2:2:C:A"]
    assert (tmp_path / "aligned.unrelated.bim").read_text(encoding="utf-8") == "keep\n"
    assert (tmp_path / "aligned.unrelated.bed").read_text(encoding="utf-8") == "keep\n"
    assert (tmp_path / "aligned.unrelated.fam").read_text(encoding="utf-8") == "keep\n"


def test_project_payload_bfile_passes_through_target_fam_and_sample_id_mode(tmp_path):
    source_prefix = tmp_path / "source"
    source_bim = source_prefix.with_suffix(".bim")
    vmap = tmp_path / "target.vmap"
    target_fam = tmp_path / "target.fam"
    output = tmp_path / "aligned"

    write_bfile(source_prefix, ["1\trs1\t0\t1\tG\tA"], ["F1 S1 0 0 1 -9", "F2 S2 0 0 2 -9"], [[0, 3]])
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])
    write_lines(target_fam, ["X9 S2 0 0 2 -9", "X8 S1 0 0 1 -9"])

    result = run_py(
        "project_payload.py",
        "--input",
        source_bim,
        "--input-format",
        "bfile",
        "--vmap",
        vmap,
        "--target-fam",
        target_fam,
        "--sample-id-mode",
        "iid",
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert "--target-fam" in result.stderr
    assert "--sample-id-mode iid" in result.stderr
    assert read_bed_matrix(output.with_suffix(".bed"), 2, 1) == [[3, 0]]


def test_project_payload_bfile_sample_axis_union_synthesizes_target_fam_from_referenced_shards_only(tmp_path):
    source_template = tmp_path / "source.@.bim"
    vmap = tmp_path / "target.vmap"
    output = tmp_path / "aligned.@"

    write_bfile(tmp_path / "source_1", ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 0 -9", "S2 S2 0 0 2 -9"], [[0, 3]])
    write_bfile(tmp_path / "source_2", ["2\trs2\t0\t200\tC\tT"], ["S2 S2 0 0 2 -9", "S3 S3 0 0 1 -9"], [[2, 0]])
    write_bfile(tmp_path / "source_3", ["3\trs3\t0\t300\tG\tT"], ["S4 S4 0 0 1 -9"], [[0]])
    for stem in ("source_1", "source_2", "source_3"):
        src = tmp_path / stem
        token = stem.split("_", 1)[1]
        for ext in (".bed", ".bim", ".fam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_vmap(
        vmap,
        ["1\t100\tt1\tA\tG\t1\t0\tidentity", "2\t200\tt2\tC\tT\t2\t0\tidentity"],
        genome_build="GRCh37",
    )

    result = run_py(
        "project_payload.py",
        "--input",
        source_template,
        "--input-format",
        "bfile",
        "--vmap",
        vmap,
        "--sample-axis",
        "union",
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    target_samples = tmp_path / "aligned.all_targets.target_samples.fam"
    assert target_samples.exists()
    assert "--target-fam" in result.stderr
    assert "S4" not in target_samples.read_text(encoding="utf-8")
    assert target_samples.read_text(encoding="utf-8") == "S1 S1 0 0 0 -9\nS2 S2 0 0 2 -9\nS3 S3 0 0 1 -9\n"
    assert (tmp_path / "aligned.1.fam").read_text(encoding="utf-8") == target_samples.read_text(encoding="utf-8")
    assert read_bed_matrix(tmp_path / "aligned.1.bed", 3, 1) == [[0, 3, 1]]
    assert read_bed_matrix(tmp_path / "aligned.2.bed", 3, 1) == [[1, 2, 0]]


def test_project_payload_bfile_passes_through_native_sample_axis_and_skip_ploidy(tmp_path):
    source_template = tmp_path / "source.@.bim"
    vmap = tmp_path / "target.vmap"
    output = tmp_path / "aligned.@"

    write_bfile(tmp_path / "source_1", ["1\trs1\t0\t100\tA\tG"], ["S1 S1 0 0 1 -9", "S2 S2 0 0 2 -9"], [[0, 3]])
    write_bfile(tmp_path / "source_X", ["X\trsx\t0\t5000000\tC\tT"], ["S1 S1 0 0 1 -9"], [[2]])
    for stem, token in (("source_1", "1"), ("source_X", "X")):
        src = tmp_path / stem
        for ext in (".bed", ".bim", ".fam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_vmap(
        vmap,
        ["1\t100\tt1\tA\tG\t1\t0\tidentity", "X\t5000000\ttx\tC\tT\tX\t0\tidentity"],
        genome_build="GRCh37",
    )

    result = run_py(
        "project_payload.py",
        "--input",
        source_template,
        "--input-format",
        "bfile",
        "--vmap",
        vmap,
        "--sample-axis",
        "native",
        "--skip-ploidy-check",
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert "--sample-axis native" in result.stderr
    assert "--skip-ploidy-check" in result.stderr
    assert (tmp_path / "aligned.1.fam").read_text(encoding="utf-8") == (tmp_path / "source.1.fam").read_text(encoding="utf-8")
    assert (tmp_path / "aligned.X.fam").read_text(encoding="utf-8") == (tmp_path / "source.X.fam").read_text(encoding="utf-8")


def test_project_payload_bfile_sample_axis_union_warns_and_noops_for_non_sharded_input(tmp_path):
    source_prefix = tmp_path / "source"
    source_bim = source_prefix.with_suffix(".bim")
    vmap = tmp_path / "target.vmap"
    output = tmp_path / "aligned"

    write_bfile(source_prefix, ["1\trs1\t0\t1\tG\tA"], ["S1 S1 0 0 1 -9"], [[0]])
    write_vmap(vmap, ["1\t1\trs1\tG\tA\t.\t0\tidentity"])

    result = run_py(
        "project_payload.py",
        "--input",
        source_bim,
        "--input-format",
        "bfile",
        "--vmap",
        vmap,
        "--sample-axis",
        "union",
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert "--sample-axis union is a no-op for non-sharded source input" in result.stderr
    assert not (tmp_path / "aligned.target_samples.fam").exists()

from __future__ import annotations

from pathlib import Path

import numpy as np
import pgenlib

from utils import run_py, write_json, write_lines


def write_psam(path: Path, sample_lines: list[str], *, header: str = "#IID") -> None:
    write_lines(path, [header, *sample_lines])


def write_pvar(path: Path, variant_lines: list[str]) -> None:
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT", *variant_lines]
    write_lines(path, lines)


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


def read_pfile_alleles_and_phase(prefix: Path, n_samples: int, variant_idx: int) -> tuple[list[int], list[int]]:
    reader = pgenlib.PgenReader(str(prefix.with_suffix(".pgen")).encode("utf-8"), n_samples)
    try:
        alleles = np.empty(n_samples * 2, dtype=np.int32)
        phasepresent = np.zeros(n_samples, dtype=np.uint8)
        reader.read_alleles_and_phasepresent(variant_idx, alleles, phasepresent)
        return alleles.astype(int).tolist(), phasepresent.astype(int).tolist()
    finally:
        reader.close()


def read_pfile_dosages(prefix: Path, n_samples: int, variant_idx: int) -> list[float]:
    reader = pgenlib.PgenReader(str(prefix.with_suffix(".pgen")).encode("utf-8"), n_samples)
    try:
        out = np.empty(n_samples, dtype=np.float32)
        reader.read_dosages(variant_idx, out)
        return out.astype(float).tolist()
    finally:
        reader.close()


def read_output_pvar(path: Path) -> list[list[str]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("##"):
            continue
        rows.append(line.lstrip("#").split("\t"))
    return rows


def write_vmap(path: Path, lines: list[str], *, genome_build: str = "GRCh37", contig_naming: str = "ncbi") -> None:
    write_lines(path, lines)
    write_json(
        path.with_name(path.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": genome_build, "contig_naming": contig_naming}},
    )


def test_apply_vmap_to_pfile_basic_mapped_row_without_swap(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG"],
        ["S1", "S2"],
        [{"channel": "hardcall", "alleles": [0, 0, 0, 1]}],
    )
    write_vmap(vmap, ["1\t100\tt1\tG\tA\t.\t0\tidentity"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode == 0, result.stderr
    assert read_output_pvar(out_prefix.with_suffix(".pvar")) == [
        ["CHROM", "POS", "ID", "REF", "ALT"],
        ["1", "100", "1:100:G:A", "A", "G"],
    ]
    assert read_pfile_genotypes(out_prefix, 2, 0) == [0, 1]
    assert out_prefix.with_suffix(".psam").read_text(encoding="utf-8") == source_prefix.with_suffix(".psam").read_text(encoding="utf-8")


def test_apply_vmap_to_pfile_retain_snp_id_uses_vmap_id(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["1\t100\trs_source\tA\tG"],
        ["S1"],
        [{"channel": "hardcall", "alleles": [0, 0]}],
    )
    write_vmap(vmap, ["1\t100\ttarget_id\tG\tA\t.\t0\tidentity"])

    result = run_py(
        "apply_vmap_to_pfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--retain-snp-id",
    )

    assert result.returncode == 0, result.stderr
    assert read_output_pvar(out_prefix.with_suffix(".pvar")) == [
        ["CHROM", "POS", "ID", "REF", "ALT"],
        ["1", "100", "target_id", "A", "G"],
    ]


def test_apply_vmap_to_pfile_swap_updates_genotypes(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG"],
        ["S1", "S2"],
        [{"channel": "hardcall", "alleles": [0, 0, 0, 1]}],
    )
    write_vmap(vmap, ["1\t100\tt1\tA\tG\t.\t0\tswap"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode == 0, result.stderr
    assert read_pfile_genotypes(out_prefix, 2, 0) == [2, 1]


def test_apply_vmap_to_pfile_flip_behaves_like_identity_and_flip_swap_like_swap(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG"],
        ["S1", "S2"],
        [{"channel": "phase", "alleles": [0, 1, 0, 1], "phasepresent": [1, 1]}],
    )
    write_vmap(
        vmap,
        [
            "1\t100\tflip_id\tG\tA\t.\t0\tflip",
            "1\t101\tflip_swap\tA\tG\t.\t0\tflip_swap",
        ],
    )

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode == 0, result.stderr
    assert read_pfile_genotypes(out_prefix, 2, 0) == [1, 1]
    assert read_pfile_genotypes(out_prefix, 2, 1) == [1, 1]
    assert read_pfile_alleles_and_phase(out_prefix, 2, 0) == ([0, 1, 0, 1], [1, 1])
    assert read_pfile_alleles_and_phase(out_prefix, 2, 1) == ([1, 0, 1, 0], [1, 1])


def test_apply_vmap_to_pfile_unmatched_row_emits_all_missing_output(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG"],
        ["S1", "S2"],
        [{"channel": "hardcall", "alleles": [0, 0, 0, 1]}],
    )
    write_vmap(vmap, ["1\t100\tt1\tG\tA\t.\t0\tidentity", "1\t200\tmissing\tT\tC\t.\t-1\tmissing"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode == 0, result.stderr
    assert read_pfile_genotypes(out_prefix, 2, 1) == [-9, -9]


def test_apply_vmap_to_pfile_only_mapped_target_drops_unmatched_rows(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG"],
        ["S1", "S2"],
        [{"channel": "hardcall", "alleles": [0, 0, 0, 1]}],
    )
    write_vmap(vmap, ["1\t100\tt1\tG\tA\t.\t0\tidentity", "1\t200\tmissing\tT\tC\t.\t-1\tmissing"])

    result = run_py(
        "apply_vmap_to_pfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--only-mapped-target",
    )

    assert result.returncode == 0, result.stderr
    assert read_output_pvar(out_prefix.with_suffix(".pvar")) == [
        ["CHROM", "POS", "ID", "REF", "ALT"],
        ["1", "100", "1:100:G:A", "A", "G"],
    ]


def test_apply_vmap_to_pfile_rejects_retained_non_biallelic_source_row(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG,C"],
        ["S1", "S2"],
        [{"channel": "hardcall", "alleles": [0, 0, 0, 1]}],
    )
    write_vmap(vmap, ["1\t100\tt1\tG\tA\t.\t0\tidentity"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode != 0
    assert "non-biallelic" in result.stderr


def test_apply_vmap_to_pfile_allows_repeated_use_of_same_mapped_source_row(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG"],
        ["S1", "S2"],
        [{"channel": "hardcall", "alleles": [0, 0, 0, 1]}],
    )
    write_vmap(vmap, ["1\t100\tt1\tG\tA\t.\t0\tidentity", "1\t101\tt1_dup\tG\tA\t.\t0\tidentity"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode == 0, result.stderr
    assert read_pfile_genotypes(out_prefix, 2, 0) == [0, 1]
    assert read_pfile_genotypes(out_prefix, 2, 1) == [0, 1]


def test_apply_vmap_to_pfile_rejects_psam_mismatch_across_shards(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(tmp_path / "source_1", ["1\t100\trs1\tA\tG"], ["S1"], [{"channel": "hardcall", "alleles": [0, 0]}])
    write_pfile(tmp_path / "source_2", ["2\t200\trs2\tC\tT"], ["S2"], [{"channel": "hardcall", "alleles": [1, 1]}])
    for stem in ("source_1", "source_2"):
        src = tmp_path / stem
        token = stem.split("_", 1)[1]
        for ext in (".pgen", ".pvar", ".psam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_vmap(vmap, ["1\t100\tt1\tG\tA\t1\t0\tidentity", "2\t200\tt2\tT\tC\t2\t0\tidentity"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_template, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode != 0
    assert ".psam shard mismatch" in result.stderr


def test_apply_vmap_to_pfile_missing_hardcalls_and_dosages_remain_missing_after_swap(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG", "1\t200\trs2\tC\tT"],
        ["S1", "S2"],
        [
            {"channel": "hardcall", "alleles": [0, 0, -9, -9]},
            {"channel": "dosage", "dosages": [0.5, -9.0]},
        ],
    )
    write_vmap(vmap, ["1\t100\thc_swap\tA\tG\t.\t0\tswap", "1\t200\tdos_swap\tT\tC\t.\t1\tswap"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode == 0, result.stderr
    assert "incoherent supported PFILE channel availability" in result.stderr
    assert read_pfile_genotypes(out_prefix, 2, 0) == [2, -9]
    assert read_pfile_dosages(out_prefix, 2, 1) == [1.5, -9.0]


def test_apply_vmap_to_pfile_reports_haploid_incompatible_hardcalls_without_rewriting(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["X\t100\trs1\tA\tG"],
        ["S1\t1"],
        [{"channel": "hardcall", "alleles": [0, 1]}],
        psam_header="#IID\tSEX",
    )
    write_vmap(vmap, ["X\t100\tt1\tG\tA\t.\t0\tidentity"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode == 0, result.stderr
    assert "found 1 incompatible diploid-het hardcalls in haploid target ploidy" in result.stderr
    assert read_pfile_genotypes(out_prefix, 1, 0) == [1]


def test_apply_vmap_to_pfile_reports_absent_hardcalls_without_rewriting(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["Y\t100\trs1\tA\tG"],
        ["S1\t2"],
        [{"channel": "hardcall", "alleles": [0, 0]}],
        psam_header="#IID\tSEX",
    )
    write_vmap(vmap, ["Y\t100\tt1\tG\tA\t.\t0\tidentity"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode == 0, result.stderr
    assert "found 1 incompatible nonmissing hardcalls in absent target regions" in result.stderr
    assert read_pfile_genotypes(out_prefix, 1, 0) == [0]


def test_apply_vmap_to_pfile_reports_absent_dosages_without_rewriting(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["Y\t100\trs1\tA\tG"],
        ["S1\t2"],
        [{"channel": "dosage", "dosages": [1.5]}],
        psam_header="#IID\tSEX",
    )
    write_vmap(vmap, ["Y\t100\tt1\tG\tA\t.\t0\tidentity"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode == 0, result.stderr
    assert "found 1 incompatible nonmissing dosages in absent target regions" in result.stderr
    assert read_pfile_dosages(out_prefix, 1, 0) == [1.5]


def test_apply_vmap_to_pfile_all_missing_row_is_coherent_for_available_channels(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG"],
        ["S1", "S2"],
        [{"channel": "dosage", "dosages": [0.5, 1.0]}],
    )
    write_vmap(vmap, ["1\t100\tt1\tG\tA\t.\t0\tidentity", "1\t200\tmissing\tT\tC\t.\t-1\tmissing"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode == 0, result.stderr
    assert read_pfile_genotypes(out_prefix, 2, 1) == [-9, -9]
    assert read_pfile_dosages(out_prefix, 2, 1) == [-9.0, -9.0]


def test_apply_vmap_to_pfile_rejects_all_unmatched_target_rows(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG"],
        ["S1", "S2"],
        [{"channel": "hardcall", "alleles": [0, 0, 0, 1]}],
    )
    write_vmap(vmap, ["1\t200\tmissing\tT\tC\t.\t-1\tmissing"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_prefix, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode != 0
    assert "will not emit an all-missing PLINK 2 payload" in result.stderr


def test_apply_vmap_to_pfile_target_psam_reorders_and_fills_missing_for_inconsistent_shards(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    target_psam = tmp_path / "target.psam"
    write_pfile(
        tmp_path / "source_1",
        ["1\t100\trs1\tA\tG"],
        ["S1", "S2"],
        [{"channel": "hardcall", "alleles": [0, 0, 1, 1]}],
    )
    write_pfile(
        tmp_path / "source_2",
        ["2\t200\trs2\tC\tT"],
        ["S2", "S3"],
        [{"channel": "hardcall", "alleles": [0, 1, 0, 0]}],
    )
    for stem in ("source_1", "source_2"):
        src = tmp_path / stem
        token = stem.split("_", 1)[1]
        for ext in (".pgen", ".pvar", ".psam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_vmap(vmap, ["1\t100\tt1\tG\tA\t1\t0\tidentity", "2\t200\tt2\tT\tC\t2\t0\tidentity"])
    write_psam(target_psam, ["S3", "S1", "S2"])

    result = run_py(
        "apply_vmap_to_pfile.py",
        "--source-prefix",
        source_template,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--target-psam",
        target_psam,
    )

    assert result.returncode == 0, result.stderr
    assert out_prefix.with_suffix(".psam").read_text(encoding="utf-8") == target_psam.read_text(encoding="utf-8")
    assert read_pfile_genotypes(out_prefix, 3, 0) == [-9, 0, 2]
    assert read_pfile_genotypes(out_prefix, 3, 1) == [0, -9, 1]
    assert "Sample-axis reconciliation summary" in result.stderr


def test_apply_vmap_to_pfile_rejects_target_psam_duplicate_subject_keys(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    target_psam = tmp_path / "target.psam"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG"],
        ["F1\tS1", "F2\tS2"],
        [{"channel": "hardcall", "alleles": [0, 0, 0, 1]}],
        psam_header="#FID\tIID",
    )
    write_vmap(vmap, ["1\t100\tt1\tG\tA\t.\t0\tidentity"])
    write_psam(target_psam, ["X1\tS1", "X2\tS1"], header="#FID\tIID")

    result = run_py(
        "apply_vmap_to_pfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--target-psam",
        target_psam,
        "--sample-id-mode",
        "iid",
    )

    assert result.returncode != 0
    assert "duplicate subject key" in result.stderr


def test_apply_vmap_to_pfile_iid_mode_ignores_fid_for_explicit_target_psam(tmp_path):
    source_prefix = tmp_path / "source"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    target_psam = tmp_path / "target.psam"
    write_pfile(
        source_prefix,
        ["1\t100\trs1\tA\tG"],
        ["F1\tS1", "F2\tS2"],
        [{"channel": "hardcall", "alleles": [0, 0, 0, 1]}],
        psam_header="#FID\tIID",
    )
    write_vmap(vmap, ["1\t100\tt1\tG\tA\t.\t0\tidentity"])
    write_psam(target_psam, ["X9\tS2", "X8\tS1"], header="#FID\tIID")

    result = run_py(
        "apply_vmap_to_pfile.py",
        "--source-prefix",
        source_prefix,
        "--vmap",
        vmap,
        "--output-prefix",
        out_prefix,
        "--target-psam",
        target_psam,
        "--sample-id-mode",
        "iid",
    )

    assert result.returncode == 0, result.stderr
    assert read_pfile_genotypes(out_prefix, 2, 0) == [1, 0]


def test_apply_vmap_to_pfile_rejects_inconsistent_fid_header_presence_under_fid_iid(tmp_path):
    source_template = tmp_path / "source.@"
    vmap = tmp_path / "map.vmap"
    out_prefix = tmp_path / "aligned"
    write_pfile(
        tmp_path / "source_1",
        ["1\t100\trs1\tA\tG"],
        ["S1"],
        [{"channel": "hardcall", "alleles": [0, 0]}],
    )
    write_pfile(
        tmp_path / "source_2",
        ["2\t200\trs2\tC\tT"],
        ["F2\tS2"],
        [{"channel": "hardcall", "alleles": [1, 1]}],
        psam_header="#FID\tIID",
    )
    for stem in ("source_1", "source_2"):
        src = tmp_path / stem
        token = stem.split("_", 1)[1]
        for ext in (".pgen", ".pvar", ".psam"):
            src.with_suffix(ext).rename(tmp_path / f"source.{token}{ext}")
    write_vmap(vmap, ["1\t100\tt1\tG\tA\t1\t0\tidentity", "2\t200\tt2\tT\tC\t2\t0\tidentity"])

    result = run_py("apply_vmap_to_pfile.py", "--source-prefix", source_template, "--vmap", vmap, "--output-prefix", out_prefix)

    assert result.returncode != 0
    assert "inconsistent .psam FID header presence" in result.stderr

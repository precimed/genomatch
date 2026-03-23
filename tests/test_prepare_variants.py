import os
import stat
from pathlib import Path

import pytest

from utils import read_tsv, run_py_with_env, write_fasta, write_json, write_lines, write_match_config


def assert_wrote(result, path: Path) -> None:
    assert f"prepare_variants.py: wrote {path}" in result.stderr


def write_fake_bcftools(path: Path) -> None:
    write_lines(
        path,
        [
            "#!/usr/bin/env python3",
            "import sys",
            "from pathlib import Path",
            "",
            "args = sys.argv[1:]",
            "input_vcf = Path(args[1])",
            "output_vcf = Path(args[args.index('-o') + 1])",
            "chain = Path(args[args.index('-c') + 1])",
            "mapping = {}",
            "for line in chain.read_text(encoding='utf-8').splitlines():",
            "    if not line.strip():",
            "        continue",
            "    src_chrom, src_pos, dst_chrom, dst_pos = line.split('\\t')[:4]",
            "    mapping[(src_chrom, src_pos)] = (dst_chrom, dst_pos)",
            "out_lines = []",
            "for line in input_vcf.read_text(encoding='utf-8').splitlines():",
            "    if line.startswith('#'):",
            "        out_lines.append(line)",
            "        continue",
            "    parts = line.split('\\t')",
            "    key = (parts[0], parts[1])",
            "    if key not in mapping:",
            "        continue",
            "    parts[0], parts[1] = mapping[key]",
            "    out_lines.append('\\t'.join(parts))",
            "output_vcf.write_text('\\n'.join(out_lines) + '\\n', encoding='utf-8')",
        ],
    )
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


def write_failing_bcftools(path: Path) -> None:
    write_lines(
        path,
        [
            "#!/usr/bin/env python3",
            "import sys",
            "sys.stderr.write('bcftools should not have been called\\n')",
            "raise SystemExit(1)",
        ],
    )
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


def write_bim(path: Path, lines: list[str]) -> None:
    write_lines(path, lines)


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


def write_sumstats_id_lookup_metadata(path: Path) -> None:
    write_lines(
        path,
        [
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
            "col_BETA: BETA",
        ],
    )


def base_env(
    tmp_path: Path,
    *,
    grch37_sequences: dict[str, str],
    grch38_sequences: dict[str, str],
    chain_lines_37_to_38: list[str] | None = None,
) -> dict[str, str]:
    grch37 = tmp_path / "hg19.fa"
    grch38 = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    write_fasta(grch37, grch37_sequences)
    write_fasta(grch38, grch38_sequences)
    write_lines(chain37to38, chain_lines_37_to_38 or ["chr1\t1\tchr1\t3", "chr1\t2\tchr1\t4", "chr2\t2\tchr2\t5"])
    write_lines(chain38to37, ["chr1\t3\tchr1\t1", "chr1\t4\tchr1\t2", "chr2\t5\tchr2\t2"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37,
        grch38_ucsc_fasta=grch38,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    return {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)}


def test_prepare_variants_output_is_stem_and_defaults_prefix_to_output(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
        grch38_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["1\trs2\t0\t2\tC\tA", "1\trs1\t0\t1\tA\tG"])

    result = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert "import_bim.py" in result.stderr
    assert "sort_variants.py" not in result.stderr
    restrict_line = next(line for line in result.stderr.splitlines() if "restrict_build_compatible.py" in line)
    assert "--sort" in restrict_line
    assert_wrote(result, tmp_path / "prepared.vmap")
    assert (tmp_path / "prepared.imported.vmap").exists()
    assert (tmp_path / "prepared.build_compatible.vmap").exists()
    assert not (tmp_path / "prepared.sorted.vmap").exists()
    assert read_tsv(tmp_path / "prepared.vmap") == read_tsv(tmp_path / "prepared.build_compatible.vmap")


def test_prepare_variants_explicit_prefix_separates_intermediates_from_final_output(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
        grch38_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "final_output"
    prefix = tmp_path / "work" / "prepared"
    write_bim(source, ["1\trs1\t0\t1\tA\tG"])

    result = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--prefix",
        prefix,
    )

    assert result.returncode == 0, result.stderr
    assert_wrote(result, tmp_path / "final_output.vmap")
    assert (tmp_path / "work" / "prepared.imported.vmap").exists()
    assert (tmp_path / "work" / "prepared.build_compatible.vmap").exists()
    assert not (tmp_path / "work" / "prepared.sorted.vmap").exists()
    assert (tmp_path / "final_output.vmap").exists()


def test_prepare_variants_runs_liftover_when_build_differs(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
        grch38_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["chr1\trs2\t0\t2\tC\tA", "chr1\trs1\t0\t1\tA\tG"])

    result = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--dst-build",
        "GRCh38",
    )

    assert result.returncode == 0, result.stderr
    assert "liftover_build.py" in result.stderr
    assert "sort_variants.py" not in result.stderr
    restrict_line = next(line for line in result.stderr.splitlines() if "restrict_build_compatible.py" in line)
    assert "--sort" not in restrict_line
    assert (tmp_path / "prepared.lifted.vmap").exists()
    assert not (tmp_path / "prepared.sorted.vmap").exists()
    assert read_tsv(tmp_path / "prepared.vmap") == read_tsv(tmp_path / "prepared.lifted.vmap")


def test_prepare_variants_optional_strand_and_contig_restriction(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
        grch38_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["chr1\trs1\t0\t1\tA\tT", "chr2\trs2\t0\t2\tC\tA"])

    result = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--dst-build",
        "GRCh37",
        "--drop-strand-ambiguous",
        "--contigs",
        "2",
    )

    assert result.returncode == 0, result.stderr
    assert "drop_strand_ambiguous.py" in result.stderr
    assert "restrict_contigs.py" in result.stderr
    assert (tmp_path / "prepared.strand_filtered.vmap").exists()
    assert (tmp_path / "prepared.contigs.vmap").exists()
    assert read_tsv(tmp_path / "prepared.vmap") == [["2", "2", "rs2", "C", "A", ".", "1", "identity"]]


def test_prepare_variants_deferred_plink_splitx_uses_interim_plink_and_retains_splitx_stage(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "A" * 10, "chr2": "A" * 10, "chrX": "A" * 200000},
        grch38_sequences={"chr1": "T" * 10, "chr2": "T" * 10, "chrX": "T" * 200000},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["chr1\trs1\t0\t1\tG\tA", "25\trsX\t0\t100000\tC\tA"])

    result = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--dst-build",
        "GRCh37",
        "--dst-contig-naming",
        "plink_splitx",
    )

    assert result.returncode == 0, result.stderr
    command_lines = [line for line in result.stderr.splitlines() if line.startswith("+ ")]
    normalize_lines = [line for line in command_lines if "normalize_contigs.py" in line]
    update_idx = next(idx for idx, line in enumerate(command_lines) if "guess_build.py" in line)
    restrict_idx = next(idx for idx, line in enumerate(command_lines) if "restrict_build_compatible.py" in line)
    assert len(normalize_lines) == 2
    assert "--to plink" in normalize_lines[0]
    assert "--to plink_splitx" in normalize_lines[1]
    first_normalize_idx = command_lines.index(normalize_lines[0])
    second_normalize_idx = command_lines.index(normalize_lines[1])
    assert first_normalize_idx < update_idx < restrict_idx < second_normalize_idx
    assert (tmp_path / "prepared.splitx.vmap").exists()
    assert read_tsv(tmp_path / "prepared.normalized.vmap") == [
        ["1", "1", "rs1", "G", "A", ".", "0", "identity"],
        ["23", "100000", "rsX", "C", "A", ".", "1", "identity"],
    ]
    assert read_tsv(tmp_path / "prepared.vmap") == read_tsv(tmp_path / "prepared.splitx.vmap")
    assert read_tsv(tmp_path / "prepared.vmap") == [
        ["1", "1", "rs1", "G", "A", ".", "0", "identity"],
        ["25", "100000", "rsX", "C", "A", ".", "1", "identity"],
    ]


def test_prepare_variants_deferred_plink_splitx_runs_after_liftover(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "A" * 10, "chr2": "A" * 10, "chrX": "A" * 30000},
        grch38_sequences={"chr1": "T" * 10, "chr2": "T" * 10, "chrX": "T" * 30000},
        chain_lines_37_to_38=["chr1\t1\tchr1\t3", "chrX\t20000\tchrX\t20000"],
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["chr1\trs1\t0\t1\tG\tA", "25\trsX\t0\t20000\tC\tA"])

    result = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--dst-build",
        "GRCh38",
        "--dst-contig-naming",
        "plink_splitx",
    )

    assert result.returncode == 0, result.stderr
    command_lines = [line for line in result.stderr.splitlines() if line.startswith("+ ")]
    liftover_idx = next(idx for idx, line in enumerate(command_lines) if "liftover_build.py" in line)
    splitx_idx = next(idx for idx, line in enumerate(command_lines) if "normalize_contigs.py" in line and "--to plink_splitx" in line)
    assert liftover_idx < splitx_idx
    assert read_tsv(tmp_path / "prepared.normalized.vmap") == [
        ["1", "1", "rs1", "G", "A", ".", "0", "identity"],
        ["23", "20000", "rsX", "C", "A", ".", "1", "identity"],
    ]
    assert read_tsv(tmp_path / "prepared.lifted.vmap") == [
        ["1", "3", "rs1", "G", "A", ".", "0", "identity"],
    ]
    assert read_tsv(tmp_path / "prepared.lifted.vmap.qc.tsv") == [
        ["source_shard", "source_index", "source_id", "status"],
        [".", "0", "rs1", "lifted"],
        [".", "1", "rsX", "ploidy_class_changed"],
    ]
    assert read_tsv(tmp_path / "prepared.splitx.vmap") == [
        ["1", "3", "rs1", "G", "A", ".", "0", "identity"],
    ]
    assert read_tsv(tmp_path / "prepared.vmap") == read_tsv(tmp_path / "prepared.splitx.vmap")


@pytest.mark.parametrize(
    "extra_args,expected_flags",
    [
        ([], ["--allow-strand-flips", "--norm-indels", "--sort"]),
        (["--no-allow-strand-flips"], ["--norm-indels", "--sort"]),
        (["--no-norm-indels"], ["--allow-strand-flips", "--sort"]),
        (["--no-allow-strand-flips", "--no-norm-indels"], ["--sort"]),
    ],
)
def test_prepare_variants_passes_restrict_build_compatible_flags_through_exactly(tmp_path, extra_args, expected_flags):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
        grch38_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["chr1\trs1\t0\t1\tA\tG"])

    result = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--dst-build",
        "GRCh37",
        *extra_args,
    )

    assert result.returncode == 0, result.stderr
    restrict_line = next(line for line in result.stderr.splitlines() if "restrict_build_compatible.py" in line)
    for flag in ("--allow-strand-flips", "--norm-indels", "--sort"):
        assert (flag in restrict_line) is (flag in expected_flags)
    assert (tmp_path / "prepared.build_compatible.vmap").exists()
    assert not (tmp_path / "prepared.sorted.vmap").exists()
    assert read_tsv(tmp_path / "prepared.vmap") == read_tsv(tmp_path / "prepared.build_compatible.vmap")


def test_prepare_variants_sumstats_requires_metadata_and_rejects_otherwise(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
        grch38_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
    )
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    write_lines(source, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "1\t1\trs1\tA\tG\t0.5"])
    write_sumstats_metadata(metadata)

    missing = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--output",
        tmp_path / "prepared",
    )
    assert missing.returncode != 0
    assert "--sumstats-metadata is required for --input-format=sumstats" in missing.stderr

    rejected = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        tmp_path / "source.bim",
        "--input-format",
        "bim",
        "--sumstats-metadata",
        metadata,
        "--output",
        tmp_path / "prepared_bim",
    )
    assert rejected.returncode != 0
    assert "--sumstats-metadata is supported only for --input-format=sumstats" in rejected.stderr

    rejected_id_vtable = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        tmp_path / "source.bim",
        "--input-format",
        "bim",
        "--id-vtable",
        tmp_path / "lookup.vtable",
        "--output",
        tmp_path / "prepared_bim_id",
    )
    assert rejected_id_vtable.returncode != 0
    assert "--id-vtable is supported only for --input-format=sumstats" in rejected_id_vtable.stderr


def test_prepare_variants_sumstats_passes_id_vtable_through_to_importer(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
        grch38_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
    )
    source = tmp_path / "study.tsv"
    metadata = tmp_path / "study.yaml"
    id_vtable = tmp_path / "lookup.vtable"
    output = tmp_path / "prepared"
    write_lines(source, ["SNP\tEA\tOA\tBETA", "rs1\tA\tG\t0.5"])
    write_sumstats_id_lookup_metadata(metadata)
    write_lines(id_vtable, ["1\t1\trs1\tC\tT"])
    write_json(
        id_vtable.with_name(id_vtable.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "sumstats",
        "--sumstats-metadata",
        metadata,
        "--id-vtable",
        id_vtable,
        "--output",
        output,
        "--dst-build",
        "GRCh37",
    )

    assert result.returncode == 0, result.stderr
    import_line = next(line for line in result.stderr.splitlines() if "import_sumstats.py" in line)
    assert "--id-vtable" in import_line
    assert str(id_vtable) in import_line
    restrict_line = next(line for line in result.stderr.splitlines() if "restrict_build_compatible.py" in line)
    assert "--sort" in restrict_line
    assert read_tsv(tmp_path / "prepared.vmap") == read_tsv(tmp_path / "prepared.build_compatible.vmap")


def test_prepare_variants_fails_when_planned_outputs_exist_without_resume_or_force(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
        grch38_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["1\trs1\t0\t1\tA\tG"])
    write_lines(tmp_path / "prepared.imported.vmap", ["1\t1\trs1\tA\tG\t.\t0\tidentity"])
    write_json(
        tmp_path / "prepared.imported.vmap.meta.json",
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
    )

    assert result.returncode != 0
    assert "planned output or intermediate files already exist" in result.stderr
    assert str(tmp_path / "prepared.imported.vmap") in result.stderr


def test_prepare_variants_rejects_resume_and_force_together(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
        grch38_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
    )
    source = tmp_path / "source.bim"
    write_bim(source, ["1\trs1\t0\t1\tA\tG"])

    result = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        tmp_path / "prepared",
        "--resume",
        "--force",
    )

    assert result.returncode != 0
    assert "not allowed with argument" in result.stderr


def test_prepare_variants_resume_skips_existing_stage_outputs_in_order(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
        grch38_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["1\trs2\t0\t2\tC\tA", "1\trs1\t0\t1\tA\tG"])

    first = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
    )
    assert first.returncode == 0, first.stderr

    for artifact in (
        tmp_path / "prepared.vmap",
        tmp_path / "prepared.vmap.meta.json",
        tmp_path / "prepared.vmap.qc.tsv",
    ):
        if artifact.exists():
            artifact.unlink()

    second = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--resume",
    )

    assert second.returncode == 0, second.stderr
    assert "skipping import_bim.py" in second.stderr
    assert "skipping restrict_build_compatible.py" in second.stderr
    assert "skipping guess_build.py" in second.stderr
    assert_wrote(second, tmp_path / "prepared.vmap")


def test_prepare_variants_resume_fails_on_hole_in_retained_stage_chain(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
        grch38_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["1\trs2\t0\t2\tC\tA", "1\trs1\t0\t1\tA\tG"])

    first = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--drop-strand-ambiguous",
    )
    assert first.returncode == 0, first.stderr

    for artifact in (
        tmp_path / "prepared.build_compatible.vmap",
        tmp_path / "prepared.build_compatible.vmap.meta.json",
        tmp_path / "prepared.vmap",
        tmp_path / "prepared.vmap.meta.json",
    ):
        artifact.unlink()

    second = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--drop-strand-ambiguous",
        "--resume",
    )

    assert second.returncode != 0
    assert "hole in retained stage outputs" in second.stderr


def test_prepare_variants_resume_noops_when_final_output_already_exists(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
        grch38_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["1\trs1\t0\t1\tA\tG"])

    first = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
    )
    assert first.returncode == 0, first.stderr

    second = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--resume",
    )

    assert second.returncode == 0, second.stderr
    assert second.stderr.strip().endswith(f"prepare_variants.py: wrote {tmp_path / 'prepared.vmap'}")
    assert "import_bim.py" not in second.stderr
    assert "restrict_build_compatible.py" not in second.stderr


def test_prepare_variants_force_deletes_outputs_first_and_reruns_cleanly(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
        grch38_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["1\trs2\t0\t2\tC\tA", "1\trs1\t0\t1\tA\tG"])

    first = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
    )
    assert first.returncode == 0, first.stderr

    write_lines(tmp_path / "prepared.vmap", ["corrupted"])

    second = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--force",
    )

    assert second.returncode == 0, second.stderr
    assert "skipping" not in second.stderr
    assert read_tsv(tmp_path / "prepared.vmap") == [
        ["1", "1", "rs1", "G", "A", ".", "1", "swap"],
        ["1", "2", "rs2", "C", "A", ".", "0", "identity"],
    ]


def test_prepare_variants_force_deletes_skipped_stage_outputs_under_prefix(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
        grch38_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["1\trs1\t0\t1\tA\tG"])
    write_lines(tmp_path / "prepared.normalized.vmap", ["stale"])
    write_json(
        tmp_path / "prepared.normalized.vmap.meta.json",
        {"object_type": "variant_map", "target": {"genome_build": "GRCh38", "contig_naming": "ncbi"}},
    )
    write_lines(tmp_path / "prepared.splitx.vmap", ["stale"])
    write_json(
        tmp_path / "prepared.splitx.vmap.meta.json",
        {"object_type": "variant_map", "target": {"genome_build": "GRCh38", "contig_naming": "plink_splitx"}},
    )

    result = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--force",
    )

    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "prepared.normalized.vmap").exists()
    assert not (tmp_path / "prepared.normalized.vmap.meta.json").exists()
    assert not (tmp_path / "prepared.splitx.vmap").exists()
    assert not (tmp_path / "prepared.splitx.vmap.meta.json").exists()


def test_prepare_variants_resume_propagates_to_liftover_reparse(tmp_path):
    env = base_env(
        tmp_path,
        grch37_sequences={"chr1": "AAAAAA", "chr2": "AAAAAA", "chrX": "AAAAAA"},
        grch38_sequences={"chr1": "TTTTTT", "chr2": "TTTTTT", "chrX": "TTTTTT"},
    )
    source = tmp_path / "source.bim"
    output = tmp_path / "prepared"
    write_bim(source, ["chr1\trs2\t0\t2\tC\tA", "chr1\trs1\t0\t1\tA\tG"])

    first = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--dst-build",
        "GRCh38",
    )
    assert first.returncode == 0, first.stderr

    for artifact in (
        tmp_path / "prepared.lifted.vmap",
        tmp_path / "prepared.lifted.vmap.meta.json",
        tmp_path / "prepared.lifted.vmap.qc.tsv",
        tmp_path / "prepared.vmap",
        tmp_path / "prepared.vmap.meta.json",
        tmp_path / "prepared.vmap.qc.tsv",
    ):
        if artifact.exists():
            artifact.unlink()
    write_failing_bcftools(tmp_path / "bcftools")

    second = run_py_with_env(
        "prepare_variants.py",
        env,
        "--input",
        source,
        "--input-format",
        "bim",
        "--output",
        output,
        "--dst-build",
        "GRCh38",
        "--resume",
    )

    assert second.returncode == 0, second.stderr
    assert "liftover_build.py --input" in second.stderr
    assert "--target-build GRCh38 --resume" in second.stderr
    assert "liftover_build.py: skipping bcftools +liftover" in second.stderr
    assert_wrote(second, tmp_path / "prepared.vmap")

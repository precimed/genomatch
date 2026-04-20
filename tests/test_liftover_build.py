import json
import os
import stat
from pathlib import Path
import pytest

from utils import read_tsv, run_py_with_env, write_fasta, write_json, write_lines, write_match_config


def write_fake_bcftools(path: Path) -> None:
    write_lines(
        path,
        [
            "#!/usr/bin/env python3",
            "import sys",
            "from pathlib import Path",
            "",
            "COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}",
            "args = sys.argv[1:]",
            "input_vcf = Path(args[1])",
            "output_vcf = Path(args[args.index('-o') + 1])",
            "chain = Path(args[args.index('-c') + 1])",
            "mapping = {}",
            "for line in chain.read_text(encoding='utf-8').splitlines():",
            "    if not line.strip():",
            "        continue",
            "    parts = line.split('\\t')",
            "    src_chrom, src_pos, dst_chrom, dst_pos = parts[:4]",
            "    extra = {item.split('=', 1)[0]: item.split('=', 1)[1] for item in parts[4:] if '=' in item}",
            "    mapping[(src_chrom, src_pos)] = (dst_chrom, dst_pos, extra)",
            "out_lines = []",
            "for line in input_vcf.read_text(encoding='utf-8').splitlines():",
            "    if line.startswith('#'):",
            "        out_lines.append(line)",
            "        continue",
            "    parts = line.split('\\t')",
            "    key = (parts[0], parts[1])",
            "    if key not in mapping:",
            "        continue",
            "    new_chrom, new_pos, extra = mapping[key]",
            "    parts[0] = new_chrom",
            "    parts[1] = new_pos",
            "    info = [] if parts[7] == '.' else parts[7].split(';')",
            "    if extra.get('flip') == '1':",
            "        parts[3] = ''.join(COMPLEMENT[base] for base in reversed(parts[3]))",
            "        parts[4] = ''.join(COMPLEMENT[base] for base in reversed(parts[4]))",
            "        info.append('FLIP')",
            "    swap = extra.get('swap')",
            "    if swap is not None:",
            "        if swap == '1':",
            "            parts[3], parts[4] = parts[4], parts[3]",
            "        info.append(f'SWAP={swap}')",
            "    parts[7] = ';'.join(info) if info else '.'",
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


def write_multiallelic_bcftools(path: Path) -> None:
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
            "out_lines = []",
            "for line in input_vcf.read_text(encoding='utf-8').splitlines():",
            "    if line.startswith('#'):",
            "        out_lines.append(line)",
            "        continue",
            "    parts = line.split('\\t')",
            "    parts[1] = '5'",
            "    parts[4] = 'G,A'",
            "    out_lines.append('\\t'.join(parts))",
            "output_vcf.write_text('\\n'.join(out_lines) + '\\n', encoding='utf-8')",
        ],
    )
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


def write_both_multibase_bcftools(path: Path) -> None:
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
            "out_lines = []",
            "for line in input_vcf.read_text(encoding='utf-8').splitlines():",
            "    if line.startswith('#'):",
            "        out_lines.append(line)",
            "        continue",
            "    parts = line.split('\\t')",
            "    parts[1] = '5'",
            "    parts[3] = 'AA'",
            "    parts[4] = 'GG'",
            "    out_lines.append('\\t'.join(parts))",
            "output_vcf.write_text('\\n'.join(out_lines) + '\\n', encoding='utf-8')",
        ],
    )
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


def test_liftover_build_vtable_emits_vtable_and_qc(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "AAAA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t5", "chr1\t2\tchr1\t3"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1", "chr1\t3\tchr1\t2"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["1\t2\trs2\tC\tA", "1\t1\trs1\tA\tG", "1\t4\trs4\tG\tA"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "3", "rs2", "C", "A"],
        ["1", "5", "rs1", "G", "A"],
    ]
    assert out.with_name(out.name + ".bcftools_input.vcf").exists()
    assert out.with_name(out.name + ".bcftools_output.vcf").exists()
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        [".", "0", "rs2", "lifted"],
        [".", "1", "rs1", "lifted"],
        [".", "2", "rs4", "unmapped"],
    ]


def test_liftover_build_resume_reuses_retained_bcftools_output(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "AAAA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t5", "chr1\t2\tchr1\t3"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1", "chr1\t3\tchr1\t2"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["1\t2\trs2\tC\tA", "1\t1\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    first = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert first.returncode == 0, first.stderr

    out.unlink()
    out.with_name(out.name + ".meta.json").unlink()
    out.with_name(out.name + ".qc.tsv").unlink()
    write_failing_bcftools(fake_bcftools)

    second = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
        "--resume",
    )
    assert second.returncode == 0, second.stderr
    assert "skipping bcftools +liftover" in second.stderr
    assert read_tsv(out) == [
        ["1", "3", "rs2", "C", "A"],
        ["1", "5", "rs1", "G", "A"],
    ]


def test_liftover_build_vmap_preserves_original_provenance(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "base.vmap"
    out = tmp_path / "lifted.vmap"

    write_fasta(grch37_ucsc, {"chr1": "AAAA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t6", "chr1\t2\tchr1\t4"])
    write_lines(chain38to37, ["chr1\t6\tchr1\t1", "chr1\t4\tchr1\t2"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["1\t1\tt1\tA\tG\tshardA\t8\tidentity", "1\t2\tt2\tT\tA\tshardA\t9\tswap"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "4", "t2", "T", "A", "shardA", "9", "swap"],
        ["1", "6", "t1", "G", "A", "shardA", "8", "swap"],
    ]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        ["shardA", "8", "t1", "lifted"],
        ["shardA", "9", "t2", "lifted"],
    ]


def test_liftover_build_drops_rows_when_lifted_ploidy_class_changes(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chrX": "A" * 60001})
    write_fasta(grch38_ucsc, {"chrX": "A" * 10000})
    write_lines(chain37to38, ["chrX\t60001\tchrX\t10000"])
    write_lines(chain38to37, ["chrX\t10000\tchrX\t60001"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["X\t60001\trs1\tG\tA"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        [".", "0", "rs1", "ploidy_class_changed"],
    ]


def test_liftover_build_uses_bcftools_flip_and_swap_tags(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vmap"
    out = tmp_path / "lifted.vmap"

    write_fasta(grch37_ucsc, {"chr1": "AAAA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAA"})
    write_lines(
        chain37to38,
        [
            "chr1\t1\tchr1\t1",
            "chr1\t2\tchr1\t2\tflip=1",
            "chr1\t3\tchr1\t3\tswap=1",
            "chr1\t4\tchr1\t4\tflip=1\tswap=1",
        ],
    )
    write_lines(chain38to37, ["chr1\t1\tchr1\t1"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(
        source,
        [
            "1\t1\trs1\tA\tG\t.\t0\tidentity",
            "1\t2\trs2\tA\tG\t.\t1\tidentity",
            "1\t3\trs3\tA\tG\t.\t2\tidentity",
            "1\t4\trs4\tA\tG\t.\t3\tidentity",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode == 0, result.stderr
    assert [row[-1] for row in read_tsv(out)] == ["swap", "flip_swap", "identity", "flip"]


def test_liftover_build_lifts_one_side_multibase_rows_with_qc(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "AA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t3", "chr1\t2\tchr1\t4"])
    write_lines(chain38to37, ["chr1\t3\tchr1\t1", "chr1\t4\tchr1\t2"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["1\t1\trs1\tA\tAT", "1\t2\trs2\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "3", "rs1", "AT", "A"],
        ["1", "4", "rs2", "G", "A"],
    ]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        [".", "0", "rs1", "lifted"],
        [".", "1", "rs2", "lifted"],
    ]


def test_liftover_build_drops_multiallelic_liftover_output_with_qc(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "A"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t5"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_multiallelic_bcftools(fake_bcftools)
    write_lines(source, ["1\t1\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        [".", "0", "rs1", "unsupported_non_snv"],
    ]


def test_liftover_build_drops_both_multibase_liftover_output_with_qc(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "A"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t5"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_both_multibase_bcftools(fake_bcftools)
    write_lines(source, ["1\t1\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        [".", "0", "rs1", "unsupported_non_snv"],
    ]


def test_liftover_build_drops_non_primary_target_contig_with_qc(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "A"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr4_GL000008v2_random\t5"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["1\t1\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        [".", "0", "rs1", "unsupported_target_contig"],
    ]


def test_liftover_build_drops_duplicate_lifted_targets_with_qc(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vmap"
    out = tmp_path / "lifted.vmap"

    write_fasta(grch37_ucsc, {"chr1": "AA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t5", "chr1\t2\tchr1\t5"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(
        source,
        [
            "1\t1\trs1\tA\tG\tshardA\t0\tidentity",
            "1\t2\trs2\tA\tG\tshardB\t0\tidentity",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        ["shardA", "0", "rs1", "duplicate_target"],
        ["shardB", "0", "rs2", "duplicate_target"],
    ]


def test_liftover_build_accepts_ucsc_named_input(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "AAAA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t5", "chr1\t2\tchr1\t3"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1", "chr1\t3\tchr1\t2"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["chr1\t2\trs2\tC\tA", "chr1\t1\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ucsc"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["chr1", "3", "rs2", "C", "A"],
        ["chr1", "5", "rs1", "G", "A"],
    ]
    meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta == {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": "ucsc"}


def test_liftover_build_fails_on_missing_chain_entry(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "AAAA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    config.write_text(
        "\n".join(
            [
                "references:",
                "  ucsc:",
                f"    GRCh37: {{fasta: {grch37_ucsc}}}",
                f"    GRCh38: {{fasta: {grch38_ucsc}}}",
                "chain:",
                "  hg38ToHg19: missing.chain",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["1\t1\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode != 0
    assert "chain.hg19tohg38" in result.stderr.lower()


def test_liftover_build_fails_on_missing_bcftools_override(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "AAAA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t5"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_lines(source, ["1\t1\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(tmp_path / "does-not-exist")},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode != 0
    assert "bcftools not found" in result.stderr.lower()


def test_liftover_build_fails_on_missing_fasta_index(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    grch37_ucsc.write_text(">chr1\nAAAA\n", encoding="utf-8")
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t5"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["1\t1\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode != 0
    assert "fasta index not found" in result.stderr.lower()


def test_liftover_build_fails_on_wrong_strand_without_flag(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "A"})
    write_fasta(grch38_ucsc, {"chr1": "A"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t1"])
    write_lines(chain38to37, ["chr1\t1\tchr1\t1"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["1\t1\trs1\tT\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode != 0
    assert "do not match reference" in result.stderr


@pytest.mark.parametrize(
    "contig_naming,chrom,pos",
    [
        ("plink", "23", "60001"),
    ],
)
def test_liftover_build_accepts_plink_family_input(tmp_path, contig_naming, chrom, pos):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / f"src.{contig_naming}.vtable"
    out = tmp_path / f"lifted.{contig_naming}.vtable"

    prefix37 = "N" * 60000
    prefix38 = "N" * 60004
    write_fasta(grch37_ucsc, {"chrX": prefix37 + "A"})
    write_fasta(grch38_ucsc, {"chrX": prefix38 + "A"})
    write_lines(chain37to38, [f"chrX\t{pos}\tchrX\t60005"])
    write_lines(chain38to37, ["chrX\t60005\tchrX\t60001"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, [f"{chrom}\t{pos}\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": contig_naming},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [[chrom, "60005", "rs1", "G", "A"]]
    meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta == {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": contig_naming}


def test_liftover_build_rejects_plink_splitx_input(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.plink_splitx.vtable"
    out = tmp_path / "lifted.plink_splitx.vtable"

    prefix37 = "N" * 60000
    prefix38 = "N" * 60004
    write_fasta(grch37_ucsc, {"chrX": prefix37 + "A"})
    write_fasta(grch38_ucsc, {"chrX": prefix38 + "A"})
    write_lines(chain37to38, ["chrX\t60001\tchrX\t60005"])
    write_lines(chain38to37, ["chrX\t60005\tchrX\t60001"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["25\t60001\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "plink_splitx"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode != 0
    assert "does not support contig_naming=plink_splitx" in result.stderr


def test_liftover_build_fails_on_unsupported_contig_label(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "AAAA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t5"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["chrUn_gl000220\t1\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode != 0
    assert "contigs inconsistent with declared contig_naming" in result.stderr.lower()


def test_liftover_build_fails_on_unsupported_swap_annotation(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "A"})
    write_fasta(grch38_ucsc, {"chr1": "A"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t1\tswap=-1"])
    write_lines(chain38to37, ["chr1\t1\tchr1\t1"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["1\t1\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode != 0
    assert "unsupported swap annotation" in result.stderr.lower()


def test_liftover_build_drops_all_unsupported_non_snv_rows_with_qc(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"
    out = tmp_path / "lifted.vtable"

    write_fasta(grch37_ucsc, {"chr1": "AAAA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t5"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["1\t1\trs1\tAT\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "liftover_build.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)},
        "--input",
        source,
        "--output",
        out,
        "--target-build",
        "GRCh38",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        [".", "0", "rs1", "unsupported_non_snv"],
    ]


def test_liftover_build_reference_access_modes_are_consistent_and_default_to_bulk(tmp_path):
    grch37_ucsc = tmp_path / "hg19.fa"
    grch38_ucsc = tmp_path / "hg38.fa"
    chain37to38 = tmp_path / "37to38.chain"
    chain38to37 = tmp_path / "38to37.chain"
    config = tmp_path / "config.yaml"
    fake_bcftools = tmp_path / "bcftools"
    source = tmp_path / "src.vtable"

    write_fasta(grch37_ucsc, {"chr1": "AAAA"})
    write_fasta(grch38_ucsc, {"chr1": "AAAAAA"})
    write_lines(chain37to38, ["chr1\t1\tchr1\t5", "chr1\t2\tchr1\t3"])
    write_lines(chain38to37, ["chr1\t5\tchr1\t1", "chr1\t3\tchr1\t2"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37_ucsc,
        grch38_ucsc_fasta=grch38_ucsc,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )
    write_fake_bcftools(fake_bcftools)
    write_lines(source, ["1\t2\trs2\tC\tA", "1\t1\trs1\tA\tG", "1\t4\trs4\tG\tA"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    def run_mode(mode: str | None):
        out = tmp_path / f"lifted_{mode or 'default'}.vtable"
        env = {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(fake_bcftools)}
        if mode is not None:
            env["MATCH_REFERENCE_ACCESS_MODE"] = mode
        result = run_py_with_env(
            "liftover_build.py",
            env,
            "--input",
            source,
            "--output",
            out,
            "--target-build",
            "GRCh38",
        )
        assert result.returncode == 0, result.stderr
        return (
            read_tsv(out),
            read_tsv(out.with_name(out.name + ".qc.tsv")),
            out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"),
        )

    default_rows, default_qc, default_meta = run_mode(None)
    bulk_rows, bulk_qc, bulk_meta = run_mode("BULK")
    legacy_rows, legacy_qc, legacy_meta = run_mode("legacy")

    assert default_rows == bulk_rows == legacy_rows
    assert default_qc == bulk_qc == legacy_qc
    assert default_meta == bulk_meta == legacy_meta

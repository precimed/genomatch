import json
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

from genomatch.vmap_restrict_build_compatible import run_bcftools_norm_check_ref_x_workaround
from utils import read_tsv, run_py_with_env, write_fasta, write_json, write_lines, write_match_config


def write_fake_norm_bcftools(path, behavior_path, log_path):
    write_lines(
        path,
        [
            "#!/usr/bin/env python3",
            "import json",
            "import sys",
            "from pathlib import Path",
            "",
            "args = sys.argv[1:]",
            f"Path({str(log_path)!r}).write_text(' '.join(args) + '\\n', encoding='utf-8')",
            "if not args or args[0] != 'norm':",
            "    sys.stderr.write('expected bcftools norm\\n')",
            "    raise SystemExit(1)",
            "check_mode = args[args.index('-c') + 1]",
            f"behavior = json.loads(Path({str(behavior_path)!r}).read_text(encoding='utf-8'))",
            "if behavior.get('__fail__'):",
            "    sys.stderr.write(behavior['__fail__'])",
            "    raise SystemExit(1)",
            "input_vcf = Path(args[-1])",
            "output_vcf = Path(args[args.index('-o') + 1])",
            "output_lines = []",
            "for line in input_vcf.read_text(encoding='utf-8').splitlines():",
            "    if line.startswith('#'):",
            "        output_lines.append(line)",
            "        continue",
            "    chrom, pos, row_id, ref, alt, qual, filt, info = line.split('\\t')[:8]",
            "    spec = behavior.get(row_id)",
            "    if spec is None:",
            "        output_lines.append(line)",
            "        continue",
            "    action = spec.get('action', 'records')",
            "    if action == 'drop':",
            "        if check_mode == 'w':",
            "            sys.stderr.write(f'REF_MISMATCH\\t{chrom}\\t{pos}\\t{ref}\\t{alt}\\n')",
            "            output_lines.append(line)",
            "        continue",
            "    for record in spec.get('records', []):",
            "        output_lines.append('\\t'.join([",
            "            record.get('chrom', chrom),",
            "            str(record.get('pos', pos)),",
            "            record.get('id', row_id),",
            "            record.get('ref', ref),",
            "            record.get('alt', alt),",
            "            '.',",
            "            'PASS',",
            "            '.',",
            "        ]))",
            "output_vcf.write_text('\\n'.join(output_lines) + '\\n', encoding='utf-8')",
        ],
    )
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


def restrict_env_with_fake_norm(tmp_path, *, behavior, grch37_sequence="AAAAAA", grch38_sequence="CCCCCC"):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    bcftools = tmp_path / "bcftools"
    behavior_path = tmp_path / "bcftools_behavior.json"
    log_path = tmp_path / "bcftools.log"
    write_fasta(grch37, {"chr1": grch37_sequence})
    write_fasta(grch38, {"chr1": grch38_sequence})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    behavior_path.write_text(json.dumps(behavior, sort_keys=True), encoding="utf-8")
    write_fake_norm_bcftools(bcftools, behavior_path, log_path)
    env = {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(bcftools)}
    return env, log_path


def norm_debug_vcf_path(output_path, *, check_mode, kind):
    return output_path.with_name(output_path.name + f".bcftools_norm_{kind}_{check_mode}.vcf")


def norm_debug_log_path(output_path, *, check_mode):
    return output_path.with_name(output_path.name + f".bcftools_norm_output_{check_mode}.log")


def write_issue_2427_fixture(tmp_path):
    fasta = tmp_path / "ref.fa"
    vcf = tmp_path / "in.vcf"
    write_fasta(fasta, {"chr1": "AAA"})
    write_lines(
        vcf,
        [
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t1\tbad\tC\tA\t.\tPASS\t.",
            "chr1\t2\tgood\tA\tAT\t.\tPASS\t.",
        ],
    )
    return fasta, vcf


def test_bcftools_issue_2427_check_ref_x_segfault_reproduced(tmp_path):
    """Sentinel: this should pass only while upstream bcftools issue #2427 exists.

    Once bcftools fixes the segfault, this test should fail and signal that the
    workaround can be removed.
    """
    bcftools = shutil.which("bcftools")
    samtools = shutil.which("samtools")
    if bcftools is None or samtools is None:
        pytest.skip("bcftools and samtools are required for the issue #2427 sentinel test")
    fasta, vcf = write_issue_2427_fixture(tmp_path)
    subprocess.run([samtools, "faidx", str(fasta)], check=True)
    output_vcf = tmp_path / "out.vcf"
    result = subprocess.run(
        [bcftools, "norm", "-f", str(fasta), "-c", "x", "-o", str(output_vcf), "-Ov", str(vcf)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode in {-11, 139}


def test_bcftools_issue_2427_workaround_with_check_ref_w_matches_expected_filtering(tmp_path):
    bcftools = shutil.which("bcftools")
    samtools = shutil.which("samtools")
    if bcftools is None or samtools is None:
        pytest.skip("bcftools and samtools are required for the issue #2427 workaround test")
    fasta, vcf = write_issue_2427_fixture(tmp_path)
    subprocess.run([samtools, "faidx", str(fasta)], check=True)
    output_vcf = tmp_path / "out.vcf"
    log_path = tmp_path / "out.log"
    run_bcftools_norm_check_ref_x_workaround(vcf, output_vcf, fasta, log_path=log_path)
    retained_records = [
        line.rstrip("\n")
        for line in output_vcf.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert retained_records == ["chr1\t2\tgood\tA\tAT\t.\tPASS\t."]
    log_text = log_path.read_text(encoding="utf-8")
    assert "[check_ref_w_pass1]" in log_text
    assert "[check_ref_w_pass2]" in log_text
    assert "REF_MISMATCH\tchr1\t1\tC\tA" in log_text


def test_restrict_build_compatible_vmap_uses_config_and_preserves_provenance(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_fasta(grch37, {"chr1": "AGT"})
    write_fasta(grch38, {"chr1": "CCC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(
        source,
        [
            "1\t1\tt1\tA\tC\tshardA\t5\tidentity",
            "1\t2\tt2\tT\tC\tshardA\t6\tidentity",
            "1\t3\tt3\tC\tG\tshardA\t7\tidentity",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
        "--allow-strand-flips",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "1", "t1", "C", "A", "shardA", "5", "swap"],
        ["1", "2", "t2", "A", "G", "shardA", "6", "flip"],
    ]
    summary = json.loads(result.stdout)
    assert summary["input"] == str(source)
    assert summary["output"] == str(out)
    assert summary["object_type"] == "variant_map"
    assert summary["genome_build"] == "GRCh37"
    assert summary["contig_naming"] == "ncbi"
    assert summary["normalization"] == "ncbi->ucsc"
    assert summary["allow_strand_flips"] is True
    assert summary["input_rows"] == 3
    assert summary["output_rows"] == 2
    assert summary["dropped_rows"] == 1
    assert summary["flip_rows"] == 1
    meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta["target"]["genome_build"] == "GRCh37"


def test_restrict_build_compatible_sort_sorts_retained_rows_into_declared_coordinate_order(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_fasta(grch37, {"chr1": "TA"})
    write_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(
        source,
        [
            "1\t2\tt2\tC\tA\tshardA\t1\tidentity",
            "1\t1\tt1\tG\tT\tshardA\t0\tidentity",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
        "--sort",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "1", "t1", "G", "T", "shardA", "0", "identity"],
        ["1", "2", "t2", "C", "A", "shardA", "1", "identity"],
    ]


def test_restrict_build_compatible_vtable_emits_vtable_and_keeps_reference_anchored_non_snv(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_fasta(grch37, {"chr1": "AG"})
    write_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tA\tAT", "1\t2\tt2\tAT\tCG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "1", "t1", "AT", "A"]]


def test_restrict_build_compatible_vtable_can_flip_reference_anchored_non_snv(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_fasta(grch37, {"chr1": "AG"})
    write_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tT\tCG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
        "--allow-strand-flips",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "1", "t1", "CG", "A"]]


def test_restrict_build_compatible_requires_known_build(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_fasta(grch37, {"chr1": "AG"})
    write_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
    )
    assert result.returncode != 0
    assert "must be known" in result.stderr.lower()


def test_restrict_build_compatible_drops_rows_where_both_alleles_are_longer_than_one_base(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_fasta(grch37, {"chr1": "AG"})
    write_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tAT\tCG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
        "--allow-strand-flips",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []


def test_restrict_build_compatible_accepts_ucsc_named_input(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_fasta(grch37, {"chr1": "AG"})
    write_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["chr1\t1\tt1\tA\tC", "chr1\t2\tt2\tT\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ucsc"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["chr1", "1", "t1", "C", "A"]]

def test_restrict_build_compatible_vmap_keeps_reference_second_rows_unchanged(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_fasta(grch37, {"chr1": "A"})
    write_fasta(grch38, {"chr1": "C"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tG\tA\tshardA\t0\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
        "--allow-strand-flips",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "1", "t1", "G", "A", "shardA", "0", "identity"]]


def test_restrict_build_compatible_vmap_swaps_when_reference_is_in_a1(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_fasta(grch37, {"chr1": "A"})
    write_fasta(grch38, {"chr1": "C"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tA\tG\tshardA\t0\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
        "--allow-strand-flips",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "1", "t1", "G", "A", "shardA", "0", "swap"]]


def test_restrict_build_compatible_vmap_can_flip_and_swap_to_reference_in_a2(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_fasta(grch37, {"chr1": "A"})
    write_fasta(grch38, {"chr1": "C"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tT\tC\tshardA\t0\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
        "--allow-strand-flips",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "1", "t1", "G", "A", "shardA", "0", "flip_swap"]]


def test_restrict_build_compatible_vmap_drops_duplicate_targets_to_qc(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_fasta(grch37, {"chr1": "A"})
    write_fasta(grch38, {"chr1": "C"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(
        source,
        [
            "1\t1\tt1\tA\tG\tshardA\t0\tidentity",
            "1\t1\tt2\tG\tA\tshardB\t1\tidentity",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        ["shardA", "0", "t1", "duplicate_target"],
        ["shardB", "1", "t2", "duplicate_target"],
    ]
    summary = json.loads(result.stdout)
    assert summary["duplicate_target_rows"] == 2


def test_restrict_build_compatible_rejects_unknown_contigs(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_fasta(grch37, {"chr1": "AG"})
    write_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["unknown\t1\tt1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
    )
    assert result.returncode != 0
    assert "normalize_contigs.py" in result.stderr


def test_restrict_build_compatible_rejects_rows_inconsistent_with_declared_contig_naming(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_fasta(grch37, {"chr1": "AG"})
    write_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["chr1\t1\tt1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
    )
    assert result.returncode != 0
    assert "normalize_contigs.py" in result.stderr


def test_restrict_build_compatible_fails_on_missing_config_section(tmp_path):
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    config.write_text("chain:\n  hg19ToHg38: missing.chain\n", encoding="utf-8")
    write_lines(source, ["1\t1\tt1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
    )
    assert result.returncode != 0
    assert "references" in result.stderr


def test_restrict_build_compatible_fails_on_missing_fasta_index(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    grch37.write_text(">chr1\nAG\n", encoding="utf-8")
    write_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
    )
    assert result.returncode != 0
    assert "fasta index not found" in result.stderr.lower()


def test_restrict_build_compatible_fails_on_bad_relative_fasta_path(tmp_path):
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    config.write_text(
        "\n".join(
            [
                "references:",
                "  ucsc:",
                "    GRCh37:",
                "      fasta: missing/GRCh37.fa",
                "    GRCh38:",
                "      fasta: missing/GRCh38.fa",
                "chain:",
                "  hg19ToHg38: missing.chain",
                "  hg38ToHg19: missing.chain",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    write_lines(source, ["1\t1\tt1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
    )
    assert result.returncode != 0
    assert "internal ucsc reference fasta" in result.stderr.lower()


def test_restrict_build_compatible_fails_on_non_acgt_alleles(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_fasta(grch37, {"chr1": "AG"})
    write_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tA\tN"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
    )
    assert result.returncode != 0
    assert "invalid allele code" in result.stderr.lower()


@pytest.mark.parametrize(
    "contig_naming,chrom,pos",
    [
        ("plink", "23", "60001"),
        ("plink_splitx", "25", "60001"),
    ],
)
def test_restrict_build_compatible_accepts_plink_family_input(tmp_path, contig_naming, chrom, pos):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / f"base.{contig_naming}.vtable"
    out = tmp_path / f"restricted.{contig_naming}.vtable"
    prefix = "N" * 60000
    write_fasta(grch37, {"chrX": prefix + "AC"})
    write_fasta(grch38, {"chrX": prefix + "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, [f"{chrom}\t{pos}\tt1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": contig_naming},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config)},
        "--source",
        source,
        "--output",
        out,
    )
    assert result.returncode == 0, result.stderr


def test_restrict_build_compatible_norm_indels_rejects_plink_splitx_input(tmp_path):
    env, _log_path = restrict_env_with_fake_norm(tmp_path, behavior={})
    source = tmp_path / "base.plink_splitx.vtable"
    out = tmp_path / "restricted.plink_splitx.vtable"
    write_lines(source, [f"25\t60001\tt1\tA\tC"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "plink_splitx"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )
    assert result.returncode != 0
    assert "does not support contig_naming=plink_splitx" in result.stderr


def test_restrict_build_compatible_norm_indels_leaves_pure_snps_on_the_non_bcftools_path(tmp_path):
    env, log_path = restrict_env_with_fake_norm(tmp_path, behavior={"__fail__": "bcftools should not have been called\n"})
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t1\tt1\tA\tG\tshardA\t0\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "1", "t1", "G", "A", "shardA", "0", "swap"]]
    assert not log_path.exists()


def test_restrict_build_compatible_norm_indels_normalizes_one_multibase_branch(tmp_path):
    env, log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={
            "branch2_row0": {
                "records": [{"chrom": "chr1", "pos": 1, "ref": "AA", "alt": "A"}],
            }
        },
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t2\tt1\tA\tAT\tshardA\t5\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )

    assert result.returncode == 0, result.stderr
    assert "-c e" in log_path.read_text(encoding="utf-8")
    assert read_tsv(out) == [["1", "1", "t1", "A", "AA", "shardA", "5", "swap"]]


def test_restrict_build_compatible_norm_indels_drops_rows_when_ploidy_class_changes(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    bcftools = tmp_path / "bcftools"
    behavior_path = tmp_path / "bcftools_behavior.json"
    log_path = tmp_path / "bcftools.log"
    write_fasta(grch37, {"chrX": "A" * 60001})
    write_fasta(grch38, {"chrX": "C" * 60001})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    behavior_path.write_text(
        json.dumps(
            {
                "branch2_row0": {
                    "records": [{"chrom": "chrX", "pos": 60000, "ref": "A", "alt": "AT"}],
                }
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    write_fake_norm_bcftools(bcftools, behavior_path, log_path)
    env = {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(bcftools)}
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["chrX\t60001\tt1\tAT\tA\tshardA\t5\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ucsc"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        ["shardA", "5", "t1", "ploidy_class_changed"],
    ]


def test_restrict_build_compatible_norm_indels_retains_branch2_vcfs_and_clears_stale_outputs(tmp_path):
    env, _log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={
            "branch2_row0": {
                "records": [{"chrom": "chr1", "pos": 1, "ref": "AA", "alt": "A"}],
            }
        },
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t2\tt1\tA\tAT\tshardA\t5\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )
    for check_mode in ("e", "x"):
        for kind in ("input", "output"):
            norm_debug_vcf_path(out, check_mode=check_mode, kind=kind).write_text("stale\n", encoding="utf-8")
        norm_debug_log_path(out, check_mode=check_mode).write_text("stale\n", encoding="utf-8")

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )

    assert result.returncode == 0, result.stderr
    branch2_input = norm_debug_vcf_path(out, check_mode="e", kind="input")
    branch2_output = norm_debug_vcf_path(out, check_mode="e", kind="output")
    assert branch2_input.exists()
    assert branch2_output.exists()
    assert norm_debug_log_path(out, check_mode="e").exists()
    assert "branch2_row0" in branch2_input.read_text(encoding="utf-8")
    assert "stale" not in branch2_input.read_text(encoding="utf-8")
    assert "stale" not in branch2_output.read_text(encoding="utf-8")
    assert "returncode: 0" in norm_debug_log_path(out, check_mode="e").read_text(encoding="utf-8")
    assert not norm_debug_vcf_path(out, check_mode="x", kind="input").exists()
    assert not norm_debug_vcf_path(out, check_mode="x", kind="output").exists()
    assert not norm_debug_log_path(out, check_mode="x").exists()


def test_restrict_build_compatible_norm_indels_composes_flip_swap_before_normalization(tmp_path):
    env, log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={
            "branch2_row0": {
                "records": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "CG"}],
            }
        },
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t1\tt1\tT\tCG\tshardA\t7\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--allow-strand-flips",
        "--norm-indels",
    )

    assert result.returncode == 0, result.stderr
    assert "-c e" in log_path.read_text(encoding="utf-8")
    assert read_tsv(out) == [["1", "1", "t1", "CG", "A", "shardA", "7", "flip_swap"]]


@pytest.mark.parametrize(
    "survivor_id,expected_allele_op",
    [
        ("branch3_row0_identity", "identity"),
        ("branch3_row0_swap", "swap"),
    ],
)
def test_restrict_build_compatible_norm_indels_branch3_two_candidate_logic(tmp_path, survivor_id, expected_allele_op):
    env, log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={
            "branch3_row0_identity": {"action": "drop"},
            "branch3_row0_swap": {"action": "drop"},
            survivor_id: {
                "records": [{"chrom": "chr1", "pos": 3, "ref": "A", "alt": "T"}],
            },
        },
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t3\tt1\tAT\tAG\tshardA\t0\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )

    assert result.returncode == 0, result.stderr
    assert "-c w" in log_path.read_text(encoding="utf-8")
    assert read_tsv(out) == [["1", "3", "t1", "T", "A", "shardA", "0", expected_allele_op]]
    assert norm_debug_vcf_path(out, check_mode="x", kind="input").exists()
    assert norm_debug_vcf_path(out, check_mode="x", kind="output").exists()
    assert norm_debug_log_path(out, check_mode="x").exists()


def test_restrict_build_compatible_norm_indels_branch3_audits_ambiguous_orientation(tmp_path):
    env, _log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={
            "branch3_row0_identity": {
                "records": [{"chrom": "chr1", "pos": 3, "ref": "A", "alt": "T"}],
            },
            "branch3_row0_swap": {
                "records": [{"chrom": "chr1", "pos": 3, "ref": "A", "alt": "G"}],
            },
        },
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t3\tt1\tAT\tAG\tshardA\t0\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        ["shardA", "0", "t1", "norm_ambiguous_orientation"],
    ]


def test_restrict_build_compatible_norm_indels_branch3_audits_ref_mismatch(tmp_path):
    env, _log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={
            "branch3_row0_identity": {"action": "drop"},
            "branch3_row0_swap": {"action": "drop"},
        },
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t3\tt1\tAT\tAG\tshardA\t0\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        ["shardA", "0", "t1", "norm_ref_mismatch"],
    ]


@pytest.mark.parametrize(
    "records,expected_status",
    [
        ([{"chrom": "chr1", "pos": 2, "ref": "A", "alt": "C,G"}], "norm_multiallelic"),
        ([{"chrom": "chr1", "pos": 2, "ref": "A", "alt": "<DEL>"}], "norm_not_atcg_alleles"),
        (
            [
                {"chrom": "chr1", "pos": 2, "ref": "A", "alt": "C"},
                {"chrom": "chr1", "pos": 3, "ref": "A", "alt": "T"},
            ],
            "norm_multiple_output_records",
        ),
        ([{"chrom": "chr1", "pos": 2, "ref": "A", "alt": "A"}], "norm_identical_ref_alt_alleles"),
        ([{"chrom": "chr1", "pos": 2, "ref": "AA", "alt": "GG"}], "norm_unsupported_complex_indel"),
        ([{"chrom": "chr1", "pos": 0, "ref": "A", "alt": "C"}], "norm_invalid_position"),
        ([{"chrom": "chrUn", "pos": 2, "ref": "A", "alt": "C"}], "unsupported_target_contig"),
    ],
)
def test_restrict_build_compatible_norm_indels_audits_normalization_specific_statuses(tmp_path, records, expected_status):
    env, _log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={"branch2_row0": {"records": records}},
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t2\tt1\tA\tAT\tshardA\t5\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        ["shardA", "5", "t1", expected_status],
    ]


def test_restrict_build_compatible_norm_indels_drops_duplicate_targets_after_normalization(tmp_path):
    env, _log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={
            "branch2_row0": {"records": [{"chrom": "chr1", "pos": 1, "ref": "AA", "alt": "A"}]},
            "branch2_row1": {"records": [{"chrom": "chr1", "pos": 1, "ref": "AA", "alt": "A"}]},
        },
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(
        source,
        [
            "1\t2\tt1\tA\tAT\tshardA\t5\tidentity",
            "1\t4\tt2\tA\tAC\tshardB\t7\tidentity",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        ["shardA", "5", "t1", "duplicate_target"],
        ["shardB", "7", "t2", "duplicate_target"],
    ]


def test_restrict_build_compatible_norm_indels_fails_when_bcftools_is_missing(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_fasta(grch37, {"chr1": "AAAAAA"})
    write_fasta(grch38, {"chr1": "CCCCCC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t2\tt1\tA\tAT"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(tmp_path / "missing-bcftools")},
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )

    assert result.returncode != 0
    assert "bcftools not found" in result.stderr.lower()


def test_restrict_build_compatible_norm_indels_fails_when_bcftools_norm_fails(tmp_path):
    env, _log_path = restrict_env_with_fake_norm(tmp_path, behavior={"__fail__": "fake bcftools failure\n"})
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_lines(source, ["1\t2\tt1\tA\tAT"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
    )

    assert result.returncode != 0
    assert "bcftools norm failed" in result.stderr.lower()
    assert str(norm_debug_log_path(out, check_mode="e")) in result.stderr
    assert "fake bcftools failure" in norm_debug_log_path(out, check_mode="e").read_text(encoding="utf-8")


def test_restrict_build_compatible_reference_access_modes_are_consistent_and_default_to_bulk(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vmap"
    write_fasta(grch37, {"chr1": "AGT"})
    write_fasta(grch38, {"chr1": "CCC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(
        source,
        [
            "1\t1\tt1\tA\tC\tshardA\t5\tidentity",
            "1\t2\tt2\tT\tC\tshardA\t6\tidentity",
            "1\t3\tt3\tC\tG\tshardA\t7\tidentity",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    def run_mode(mode: str | None):
        out = tmp_path / f"restricted_{mode or 'default'}.vmap"
        env = {"MATCH_CONFIG": str(config)}
        if mode is not None:
            env["MATCH_REFERENCE_ACCESS_MODE"] = mode
        result = run_py_with_env(
            "restrict_build_compatible.py",
            env,
            "--source",
            source,
            "--output",
            out,
            "--allow-strand-flips",
        )
        assert result.returncode == 0, result.stderr
        return (
            read_tsv(out),
            json.loads(result.stdout),
            out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"),
        )

    default_rows, default_summary, default_meta = run_mode(None)
    bulk_rows, bulk_summary, bulk_meta = run_mode("BULK")
    legacy_rows, legacy_summary, legacy_meta = run_mode("legacy")

    assert default_rows == bulk_rows == legacy_rows
    assert default_meta == bulk_meta == legacy_meta
    for summary in (default_summary, bulk_summary, legacy_summary):
        summary.pop("output", None)
    assert default_summary == bulk_summary == legacy_summary

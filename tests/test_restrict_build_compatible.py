import json
import os
import stat
from pathlib import Path

import pytest

from utils import read_tsv, run_py_with_env, write_object_metadata, write_lines, write_match_config, write_primary_ucsc_fasta


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


def restrict_env_with_fake_norm(
    tmp_path,
    *,
    behavior,
    grch37_sequence="AAAAAA",
    grch38_sequence="CCCCCC",
    t2t_sequence=None,
):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    t2t = tmp_path / "T2T-CHM13v2.0.ucsc.fa"
    config = tmp_path / "config.yaml"
    bcftools = tmp_path / "bcftools"
    behavior_path = tmp_path / "bcftools_behavior.json"
    log_path = tmp_path / "bcftools.log"
    write_primary_ucsc_fasta(grch37, {"chr1": grch37_sequence})
    write_primary_ucsc_fasta(grch38, {"chr1": grch38_sequence})
    if t2t_sequence is not None:
        write_primary_ucsc_fasta(t2t, {"chr1": t2t_sequence})
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37,
        grch38_ucsc_fasta=grch38,
        t2t_chm13v2_ucsc_fasta=t2t if t2t_sequence is not None else None,
    )
    behavior_path.write_text(json.dumps(behavior, sort_keys=True), encoding="utf-8")
    write_fake_norm_bcftools(bcftools, behavior_path, log_path)
    env = {"MATCH_CONFIG": str(config), "MATCH_BCFTOOLS": str(bcftools)}
    return env, log_path


def norm_debug_vcf_path(output_path, *, check_mode, kind):
    return output_path.with_name(output_path.name + (f".bcftools_norm_batch0.{kind}.vcf" if check_mode == "e" else f".bcftools_norm_{kind}_{check_mode}.vcf"))


def norm_debug_log_path(output_path, *, check_mode):
    return output_path.with_name(output_path.name + (f".bcftools_norm_batch0.log" if check_mode == "e" else f".bcftools_norm_output_{check_mode}.log"))


def test_restrict_build_compatible_vmap_uses_config_and_preserves_provenance(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_primary_ucsc_fasta(grch37, {"chr1": "AGT"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CCC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(
        source,
        [
            "1\t1\tt1\tA\tC\tshardA\t5\tidentity",
            "1\t2\tt2\tT\tC\tshardA\t6\tidentity",
            "1\t3\tt3\tC\tG\tshardA\t7\tidentity",
        ],
    )
    write_object_metadata(
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


def test_restrict_build_compatible_accepts_t2t_build(tmp_path):
    env, _log_path = restrict_env_with_fake_norm(tmp_path, behavior={}, t2t_sequence="A")
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t1\tt1\tG\tA\tshardA\t0\tidentity"])
    write_object_metadata(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "T2T-CHM13v2.0", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "1", "t1", "G", "A", "shardA", "0", "identity"]]
    summary = json.loads(result.stdout)
    assert summary["genome_build"] == "T2T-CHM13v2.0"
    meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta["target"]["genome_build"] == "T2T-CHM13v2.0"


def test_restrict_build_compatible_t2t_fails_when_fasta_is_absent(tmp_path):
    env, _log_path = restrict_env_with_fake_norm(tmp_path, behavior={})
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t1\tt1\tG\tA\tshardA\t0\tidentity"])
    write_object_metadata(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "T2T-CHM13v2.0", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
    )

    assert result.returncode != 0
    assert "config has no builds.T2T-CHM13v2.0.ucsc_fasta entry" in result.stderr


def test_restrict_build_compatible_preserves_retained_input_order(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_primary_ucsc_fasta(grch37, {"chr1": "TA"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(
        source,
        [
            "1\t2\tt2\tC\tA\tshardA\t1\tidentity",
            "1\t1\tt1\tG\tT\tshardA\t0\tidentity",
        ],
    )
    write_object_metadata(
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
    assert read_tsv(out) == [
        ["1", "2", "t2", "C", "A", "shardA", "1", "identity"],
        ["1", "1", "t1", "G", "T", "shardA", "0", "identity"],
    ]


def test_restrict_build_compatible_vtable_emits_vtable_and_keeps_reference_anchored_non_snv(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_primary_ucsc_fasta(grch37, {"chr1": "AG"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tA\tAT", "1\t2\tt2\tAT\tCG"])
    write_object_metadata(
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
        "--include-indels",
        "--allow-indel-swaps",
    )
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "1", "t1", "AT", "A"]]


def test_restrict_build_compatible_vtable_does_not_flip_multibase_input(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_primary_ucsc_fasta(grch37, {"chr1": "AG"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tT\tCG"])
    write_object_metadata(
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
    assert read_tsv(out.with_name(out.name + ".qc.tsv"))[-1][-1] == "non_snv_input"


def test_restrict_build_compatible_requires_known_build(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_primary_ucsc_fasta(grch37, {"chr1": "AG"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tA\tC"])
    write_object_metadata(
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
    write_primary_ucsc_fasta(grch37, {"chr1": "AG"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tAT\tCG"])
    write_object_metadata(
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
    write_primary_ucsc_fasta(grch37, {"chr1": "AG"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["chr1\t1\tt1\tA\tC", "chr1\t2\tt2\tT\tC"])
    write_object_metadata(
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
    write_primary_ucsc_fasta(grch37, {"chr1": "A"})
    write_primary_ucsc_fasta(grch38, {"chr1": "C"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tG\tA\tshardA\t0\tidentity"])
    write_object_metadata(
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
    write_primary_ucsc_fasta(grch37, {"chr1": "A"})
    write_primary_ucsc_fasta(grch38, {"chr1": "C"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tA\tG\tshardA\t0\tidentity"])
    write_object_metadata(
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
    write_primary_ucsc_fasta(grch37, {"chr1": "A"})
    write_primary_ucsc_fasta(grch38, {"chr1": "C"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tT\tC\tshardA\t0\tidentity"])
    write_object_metadata(
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
    write_primary_ucsc_fasta(grch37, {"chr1": "A"})
    write_primary_ucsc_fasta(grch38, {"chr1": "C"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(
        source,
        [
            "1\t1\tt1\tA\tG\tshardA\t0\tidentity",
            "1\t1\tt2\tG\tA\tshardB\t1\tidentity",
        ],
    )
    write_object_metadata(
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
    write_primary_ucsc_fasta(grch37, {"chr1": "AG"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["unknown\t1\tt1\tA\tC"])
    write_object_metadata(
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
    write_primary_ucsc_fasta(grch37, {"chr1": "AG"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["chr1\t1\tt1\tA\tC"])
    write_object_metadata(
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
    write_object_metadata(
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
    assert "legacy match config shape" in result.stderr


def test_restrict_build_compatible_fails_on_missing_fasta_index(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    grch37.write_text(">chr1\nAG\n", encoding="utf-8")
    write_primary_ucsc_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tA\tC"])
    write_object_metadata(
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
                "builds:",
                "  GRCh37:",
                "    ucsc_fasta: missing/GRCh37.fa",
                "  GRCh38:",
                "    ucsc_fasta: missing/GRCh38.fa",
                "liftover: []",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    write_lines(source, ["1\t1\tt1\tA\tC"])
    write_object_metadata(
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
    assert "builds.grch37.ucsc_fasta" in result.stderr.lower()


def test_restrict_build_compatible_fails_on_non_acgt_alleles(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_primary_ucsc_fasta(grch37, {"chr1": "AG"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t1\tt1\tA\tN"])
    write_object_metadata(
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
    write_primary_ucsc_fasta(grch37, {"chrX": prefix + "AC"})
    write_primary_ucsc_fasta(grch38, {"chrX": prefix + "TG"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, [f"{chrom}\t{pos}\tt1\tA\tC"])
    write_object_metadata(
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
    write_object_metadata(
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
        "--include-indels",
        "--allow-indel-swaps",
    )
    assert result.returncode != 0
    assert "does not support contig_naming=plink_splitx" in result.stderr


def test_restrict_build_compatible_norm_indels_leaves_pure_snps_on_the_non_bcftools_path(tmp_path):
    env, log_path = restrict_env_with_fake_norm(tmp_path, behavior={"__fail__": "bcftools should not have been called\n"})
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t1\tt1\tA\tG\tshardA\t0\tidentity"])
    write_object_metadata(
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
        "--include-indels",
        "--allow-indel-swaps",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "1", "t1", "G", "A", "shardA", "0", "swap"]]
    assert not log_path.exists()


def test_restrict_build_compatible_norm_indels_normalizes_one_multibase_branch(tmp_path):
    env, log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={
            "row0": {
                "records": [{"chrom": "chr1", "pos": 1, "ref": "AA", "alt": "A"}],
            }
        },
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t2\tt1\tA\tAT\tshardA\t5\tidentity"])
    write_object_metadata(
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
        "--include-indels",
        "--allow-indel-swaps",
    )

    assert result.returncode == 0, result.stderr
    assert "-c e" in log_path.read_text(encoding="utf-8")
    assert read_tsv(out) == [["1", "1", "t1", "A", "AA", "shardA", "5", "swap"]]


def test_restrict_build_compatible_norm_indels_accepts_t2t_build(tmp_path):
    env, log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={
            "row0": {
                "records": [{"chrom": "chr1", "pos": 1, "ref": "AA", "alt": "A"}],
            }
        },
        t2t_sequence="AAAAAA",
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t2\tt1\tA\tAT\tshardA\t5\tidentity"])
    write_object_metadata(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "T2T-CHM13v2.0", "contig_naming": "ncbi"}},
    )

    result = run_py_with_env(
        "restrict_build_compatible.py",
        env,
        "--source",
        source,
        "--output",
        out,
        "--norm-indels",
        "--include-indels",
        "--allow-indel-swaps",
    )

    assert result.returncode == 0, result.stderr
    assert "-c e" in log_path.read_text(encoding="utf-8")
    assert read_tsv(out) == [["1", "1", "t1", "A", "AA", "shardA", "5", "swap"]]
    summary = json.loads(result.stdout)
    assert summary["genome_build"] == "T2T-CHM13v2.0"


def test_restrict_build_compatible_norm_indels_drops_rows_when_ploidy_class_changes(tmp_path):
    grch37 = tmp_path / "GRCh37.ucsc.fa"
    grch38 = tmp_path / "GRCh38.ucsc.fa"
    config = tmp_path / "config.yaml"
    bcftools = tmp_path / "bcftools"
    behavior_path = tmp_path / "bcftools_behavior.json"
    log_path = tmp_path / "bcftools.log"
    write_primary_ucsc_fasta(grch37, {"chrX": "A" * 60001})
    write_primary_ucsc_fasta(grch38, {"chrX": "C" * 60001})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    behavior_path.write_text(
        json.dumps(
            {
                "row0": {
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
    write_object_metadata(
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
        "--include-indels",
        "--allow-indel-swaps",
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
            "row0": {
                "records": [{"chrom": "chr1", "pos": 1, "ref": "AA", "alt": "A"}],
            }
        },
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t2\tt1\tA\tAT\tshardA\t5\tidentity"])
    write_object_metadata(
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
        "--include-indels",
        "--allow-indel-swaps",
    )

    assert result.returncode == 0, result.stderr
    branch2_input = norm_debug_vcf_path(out, check_mode="e", kind="input")
    branch2_output = norm_debug_vcf_path(out, check_mode="e", kind="output")
    assert branch2_input.exists()
    assert branch2_output.exists()
    assert norm_debug_log_path(out, check_mode="e").exists()
    assert "row0" in branch2_input.read_text(encoding="utf-8")
    assert "stale" not in branch2_input.read_text(encoding="utf-8")
    assert "stale" not in branch2_output.read_text(encoding="utf-8")
    assert "returncode: 0" in norm_debug_log_path(out, check_mode="e").read_text(encoding="utf-8")
    assert not norm_debug_vcf_path(out, check_mode="x", kind="input").exists()
    assert not norm_debug_vcf_path(out, check_mode="x", kind="output").exists()
    assert not norm_debug_log_path(out, check_mode="x").exists()


def test_restrict_build_compatible_norm_indels_does_not_flip_multibase_before_normalization(tmp_path):
    env, log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={
            "row0": {
                "records": [{"chrom": "chr1", "pos": 1, "ref": "A", "alt": "CG"}],
            }
        },
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t1\tt1\tT\tCG\tshardA\t7\tidentity"])
    write_object_metadata(
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
        "--include-indels",
        "--allow-indel-swaps",
    )

    assert result.returncode == 0, result.stderr
    assert not log_path.exists()
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv"))[-1][-1] == "no_reference_match"


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
        ([{"chrom": "chr1", "pos": 0, "ref": "A", "alt": "C"}], "norm_invalid_position"),
        ([{"chrom": "chrUn", "pos": 2, "ref": "A", "alt": "C"}], "unsupported_target_contig"),
    ],
)
def test_restrict_build_compatible_norm_indels_audits_normalization_specific_statuses(tmp_path, records, expected_status):
    env, _log_path = restrict_env_with_fake_norm(
        tmp_path,
        behavior={"row0": {"records": records}},
    )
    source = tmp_path / "base.vmap"
    out = tmp_path / "restricted.vmap"
    write_lines(source, ["1\t2\tt1\tA\tAT\tshardA\t5\tidentity"])
    write_object_metadata(
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
        "--include-indels",
        "--allow-indel-swaps",
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
            "row0": {"records": [{"chrom": "chr1", "pos": 1, "ref": "AA", "alt": "A"}]},
            "row1": {"records": [{"chrom": "chr1", "pos": 1, "ref": "AA", "alt": "A"}]},
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
    write_object_metadata(
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
        "--include-indels",
        "--allow-indel-swaps",
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
    write_primary_ucsc_fasta(grch37, {"chr1": "AAAAAA"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CCCCCC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(source, ["1\t2\tt1\tA\tAT"])
    write_object_metadata(
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
        "--include-indels",
        "--allow-indel-swaps",
    )

    assert result.returncode != 0
    assert "bcftools not found" in result.stderr.lower()


def test_restrict_build_compatible_norm_indels_fails_when_bcftools_norm_fails(tmp_path):
    env, _log_path = restrict_env_with_fake_norm(tmp_path, behavior={"__fail__": "fake bcftools failure\n"})
    source = tmp_path / "base.vtable"
    out = tmp_path / "restricted.vtable"
    write_lines(source, ["1\t2\tt1\tA\tAT"])
    write_object_metadata(
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
        "--include-indels",
        "--allow-indel-swaps",
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
    write_primary_ucsc_fasta(grch37, {"chr1": "AGT"})
    write_primary_ucsc_fasta(grch38, {"chr1": "CCC"})
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    write_lines(
        source,
        [
            "1\t1\tt1\tA\tC\tshardA\t5\tidentity",
            "1\t2\tt2\tT\tC\tshardA\t6\tidentity",
            "1\t3\tt3\tC\tG\tshardA\t7\tidentity",
        ],
    )
    write_object_metadata(
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


@pytest.mark.parametrize("flag", ["--sort", "--drop-duplicates"])
@pytest.mark.parametrize("suffix", ["vmap", "vtable"])
def test_restrict_build_compatible_rejects_removed_flags(tmp_path, flag, suffix):
    result = run_py_with_env(
        "restrict_build_compatible.py", {},
        "--source", tmp_path / f"source.{suffix}",
        "--output", tmp_path / f"output.{suffix}", flag,
    )
    assert result.returncode != 0
    assert f"unrecognized arguments: {flag}" in result.stderr


def test_imported_pvar_with_one_reversed_row_is_oriented_per_row(tmp_path):
    pvar = tmp_path / 'mixed.pvar'
    imported = tmp_path / 'imported.vmap'
    output = tmp_path / 'compatible.vmap'
    fasta = tmp_path / 'reference.fa'
    config = tmp_path / 'config.yaml'
    write_primary_ucsc_fasta(fasta, {'chr1': 'ACG'})
    write_match_config(config, grch37_ucsc_fasta=fasta, grch38_ucsc_fasta=fasta)
    env = {'MATCH_CONFIG': str(config)}
    write_lines(pvar, [
        '##fileformat=VCFv4.2',
        '#CHROM\tPOS\tID\tREF\tALT',
        '1\t1\tnormal1\tA\tG',
        '1\t2\treversed\tT\tC',  # FASTA is C; only this row requires swap.
        '1\t3\tnormal2\tG\tA',
    ])
    result = run_py_with_env('import_pvar.py', env, '--input', pvar,
                             '--output', imported, '--genome-build', 'GRCh37')
    assert result.returncode == 0, result.stderr
    assert read_tsv(imported) == [
        ['1', '1', 'normal1', 'G', 'A', '.', '0', 'identity'],
        ['1', '2', 'reversed', 'C', 'T', '.', '1', 'identity'],
        ['1', '3', 'normal2', 'A', 'G', '.', '2', 'identity'],
    ]
    result = run_py_with_env('restrict_build_compatible.py', env, '--source', imported, '--output', output)
    assert result.returncode == 0, result.stderr
    assert read_tsv(output) == [
        ['1', '1', 'normal1', 'G', 'A', '.', '0', 'identity'],
        ['1', '2', 'reversed', 'T', 'C', '.', '1', 'swap'],
        ['1', '3', 'normal2', 'A', 'G', '.', '2', 'identity'],
    ]
    assert not output.with_name(output.name + '.qc.tsv').exists()

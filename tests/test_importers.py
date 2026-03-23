import json

from utils import read_tsv, run_py, write_lines


def write_vtable_with_meta(path, lines, *, genome_build="GRCh37", contig_naming="ucsc"):
    write_lines(path, lines)
    path.with_name(path.name + ".meta.json").write_text(
        json.dumps(
            {
                "object_type": "variant_table",
                "genome_build": genome_build,
                **({"contig_naming": contig_naming} if contig_naming is not None else {}),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def test_import_bim_writes_vmap_and_single_file_provenance(tmp_path):
    source = tmp_path / "input.bim"
    out = tmp_path / "out.vmap"
    write_lines(source, ["chr1\trs1\t0\t100\tA\tG", "chrX\trs2\t0\t200\tC\tT"])

    result = run_py("import_bim.py", "--input", source, "--output", out, "--genome-build", "GRCh37")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["chr1", "100", "rs1", "A", "G", ".", "0", "identity"],
        ["chrX", "200", "rs2", "C", "T", ".", "1", "identity"],
    ]
    meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta["object_type"] == "variant_map"
    assert meta["target"]["genome_build"] == "GRCh37"
    assert meta["target"]["contig_naming"] == "ucsc"


def test_import_bim_creates_missing_output_directory(tmp_path):
    source = tmp_path / "input.bim"
    out = tmp_path / "nested" / "out" / "cohort.vmap"
    write_lines(source, ["chr1\trs1\t0\t100\tA\tG"])

    result = run_py("import_bim.py", "--input", source, "--output", out, "--genome-build", "GRCh37")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["chr1", "100", "rs1", "A", "G", ".", "0", "identity"]]
    meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta["target"]["genome_build"] == "GRCh37"


def test_import_bim_numeric_autosomes_infer_ncbi_without_warning(tmp_path):
    source = tmp_path / "input.bim"
    out = tmp_path / "out.vmap"
    write_lines(source, ["1\trs1\t0\t100\tA\tG", "22\trs2\t0\t200\tC\tT"])

    result = run_py("import_bim.py", "--input", source, "--output", out, "--genome-build", "GRCh37")
    assert result.returncode == 0, result.stderr
    assert "Warning:" not in result.stderr
    meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta["target"]["contig_naming"] == "ncbi"


def test_import_bim_preserves_mixed_contigs_and_omits_contig_naming(tmp_path):
    source = tmp_path / "input.bim"
    out = tmp_path / "out.vmap"
    write_lines(source, ["chr1\trs1\t0\t100\tA\tG", "2\trs2\t0\t200\tC\tT", "chrUn\trs3\t0\t300\tG\tA"])

    result = run_py("import_bim.py", "--input", source, "--output", out, "--genome-build", "GRCh37")
    assert result.returncode == 0, result.stderr
    assert "normalize_contigs.py" in result.stderr
    assert read_tsv(out) == [
        ["chr1", "100", "rs1", "A", "G", ".", "0", "identity"],
        ["2", "200", "rs2", "C", "T", ".", "1", "identity"],
        ["chrUn", "300", "rs3", "G", "A", ".", "2", "identity"],
    ]
    meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert "contig_naming" not in meta["target"]


def test_import_bim_discovers_shards_preserves_exact_source_shard_and_writes_qc(tmp_path):
    out = tmp_path / "out.vmap"
    write_lines(tmp_path / "input.chr1.bim", ["chr1\trs1\t0\t100\tA\tG", "chr2\trs2\t0\t200\tC\tT"])
    write_lines(tmp_path / "input.2.bim", ["2\trs3\t0\t300\tT\tC"])

    result = run_py("import_bim.py", "--input", tmp_path / "input.@.bim", "--output", out, "--contigs", "2")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["2", "300", "rs3", "T", "C", "2", "0", "identity"],
        ["chr2", "200", "rs2", "C", "T", "chr1", "1", "identity"],
    ]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "reason"],
        ["chr1", "0", "filtered_by_chr2use"],
    ]


def test_import_bim_discovers_xy_shard_alias_for_chr_template(tmp_path):
    out = tmp_path / "out.vmap"
    write_lines(tmp_path / "input.chrX.bim", ["23\trs_x\t0\t100\tA\tG"])
    write_lines(tmp_path / "input.chrXY.bim", ["25\trs_xy\t0\t200\tC\tT"])

    result = run_py("import_bim.py", "--input", tmp_path / "input.chr@.bim", "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["23", "100", "rs_x", "A", "G", "X", "0", "identity"],
        ["25", "200", "rs_xy", "C", "T", "XY", "0", "identity"],
    ]
    meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta["target"]["contig_naming"] == "plink_splitx"


def test_import_bim_writes_qc_for_malformed_and_non_actg_rows(tmp_path):
    source = tmp_path / "input.bim"
    out = tmp_path / "out.vmap"
    write_lines(
        source,
        [
            "1\trs_bad\t0\t100\tA",
            "1\trs_nonactg\t0\t200\tA\tN",
            "1\trs_ok\t0\t300\tC\tT",
        ],
    )

    result = run_py("import_bim.py", "--input", source, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "300", "rs_ok", "C", "T", ".", "2", "identity"]]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "reason"],
        [".", "0", "malformed_row"],
        [".", "1", "non_actg_allele"],
    ]


def test_import_vcf_drops_multiallelic_rows_to_qc(tmp_path):
    source = tmp_path / "input.vcf"
    out = tmp_path / "out.vmap"
    write_lines(
        source,
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\trs1\tG\tA,C\t.\t.\t.",
            "1\t200\trs2\tC\tA\t.\t.\t.",
        ],
    )
    result = run_py("import_vcf.py", "--input", source, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "200", "rs2", "A", "C", ".", "1", "identity"]]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "reason"],
        [".", "0", "multiallelic"],
    ]


def test_import_vcf_writes_qc_for_malformed_and_non_actg_rows(tmp_path):
    source = tmp_path / "input.vcf"
    out = tmp_path / "out.vmap"
    write_lines(
        source,
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\trs_bad",
            "1\t200\trs_nonactg\tC\tN\t.\t.\t.",
            "1\t300\trs_ok\tG\tA\t.\t.\t.",
        ],
    )

    result = run_py("import_vcf.py", "--input", source, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "300", "rs_ok", "A", "G", ".", "2", "identity"]]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "reason"],
        [".", "0", "malformed_row"],
        [".", "1", "non_actg_allele"],
    ]


def test_import_vcf_discovers_shards_preserves_exact_source_shard_and_writes_qc(tmp_path):
    out = tmp_path / "out.vmap"
    write_lines(
        tmp_path / "input.chr1.vcf",
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t100\trs1\tG\tA\t.\t.\t.",
            "chr2\t200\trs2\tC\tA\t.\t.\t.",
        ],
    )
    write_lines(
        tmp_path / "input.2.vcf",
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "2\t300\trs3\tT\tC\t.\t.\t.",
        ],
    )

    result = run_py("import_vcf.py", "--input", tmp_path / "input.@.vcf", "--output", out, "--contigs", "2")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["2", "300", "rs3", "C", "T", "2", "0", "identity"],
        ["chr2", "200", "rs2", "A", "C", "chr1", "1", "identity"],
    ]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "reason"],
        ["chr1", "0", "filtered_by_chr2use"],
    ]


def test_import_pvar_drops_multiallelic_rows_to_qc(tmp_path):
    source = tmp_path / "input.pvar"
    out = tmp_path / "out.vmap"
    write_lines(
        source,
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT",
            "1\t100\trs1\tG\tA,C",
            "1\t200\trs2\tC\tA",
        ],
    )

    result = run_py("import_pvar.py", "--input", source, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "200", "rs2", "A", "C", ".", "1", "identity"]]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "reason"],
        [".", "0", "multiallelic"],
    ]


def test_import_pvar_discovers_shards_preserves_exact_source_shard_and_writes_qc(tmp_path):
    out = tmp_path / "out.vmap"
    write_lines(
        tmp_path / "input.chr1.pvar",
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT",
            "chr1\t100\trs1\tG\tA",
            "chr2\t200\trs2\tC\tA",
        ],
    )
    write_lines(
        tmp_path / "input.2.pvar",
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT",
            "2\t300\trs3\tT\tC",
        ],
    )

    result = run_py("import_pvar.py", "--input", tmp_path / "input.@.pvar", "--output", out, "--contigs", "2")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["2", "300", "rs3", "C", "T", "2", "0", "identity"],
        ["chr2", "200", "rs2", "A", "C", "chr1", "1", "identity"],
    ]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "reason"],
        ["chr1", "0", "filtered_by_chr2use"],
    ]


def test_import_sumstats_uses_metadata_contract_and_single_file_provenance(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    out = tmp_path / "ss.vmap"
    write_lines(sumstats, ["CHR\tPOS\tSNP\tEA\tOA\tBETA", "X\t100\trs1\tA\tG\t0.2"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_SNP: SNP", "col_EffectAllele: EA", "col_OtherAllele: OA"])
    result = run_py("import_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["X", "100", "rs1", "A", "G", ".", "0", "identity"]]
    meta_payload = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta_payload["target"]["contig_naming"] == "ncbi"


def test_import_sumstats_rejects_template_input(tmp_path):
    meta = tmp_path / "ss.yaml"
    out = tmp_path / "ss.vmap"
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA"])

    result = run_py("import_sumstats.py", "--input", tmp_path / "ss.@.tsv", "--sumstats-metadata", meta, "--output", out)
    assert result.returncode != 0
    assert "does not accept '@' paths" in result.stderr


def test_import_sumstats_extracts_joined_variant_fields(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    out = tmp_path / "ss.vmap"
    write_lines(sumstats, ["VARIANT\tSNP", "1_100:A:G\trs1"])
    write_lines(
        meta,
        [
            "col_CHR: VARIANT",
            "col_POS: VARIANT",
            "col_SNP: SNP",
            "col_EffectAllele: VARIANT",
            "col_OtherAllele: VARIANT",
        ],
    )

    result = run_py("import_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "rs1", "A", "G", ".", "0", "identity"]]


def test_import_sumstats_rejects_ambiguous_metadata_keys(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    out = tmp_path / "ss.vmap"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA", "1\t100\tA\tG"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "COL_chr: CHR",
            "col_POS: POS",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
        ],
    )

    result = run_py("import_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--output", out)
    assert result.returncode != 0
    assert "ambiguous top-level keys differing only by case" in result.stderr


def test_import_sumstats_does_not_resolve_nested_column_mappings(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    out = tmp_path / "ss.vmap"
    write_lines(sumstats, ["CHR\tPOS\tEA\tOA", "1\t100\tA\tG"])
    write_lines(
        meta,
        [
            "columns:",
            "  col_CHR: CHR",
            "  col_POS: POS",
            "  col_EffectAllele: EA",
            "  col_OtherAllele: OA",
        ],
    )

    result = run_py("import_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--output", out)
    assert result.returncode != 0
    assert "missing required column mapping: col_CHR" in result.stderr


def test_import_sumstats_fails_on_headerless_input(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    out = tmp_path / "ss.vmap"
    write_lines(sumstats, ["1\t100\tA\tG"])
    write_lines(meta, ["col_CHR: CHR", "col_POS: POS", "col_EffectAllele: EA", "col_OtherAllele: OA"])

    result = run_py("import_sumstats.py", "--input", sumstats, "--sumstats-metadata", meta, "--output", out)
    assert result.returncode != 0
    assert "column not found for col_CHR" in result.stderr


def test_import_sumstats_id_vtable_enriches_coordinates_and_inherits_metadata(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    id_vtable = tmp_path / "lookup.vtable"
    out = tmp_path / "ss.vmap"
    write_lines(sumstats, ["SNP\tEA\tOA\tBETA", "rs1\tA\tG\t0.2"])
    write_lines(meta, ["col_SNP: SNP", "col_EffectAllele: EA", "col_OtherAllele: OA"])
    write_vtable_with_meta(id_vtable, ["chr1\t100\trs1\tC\tT"], genome_build="GRCh38", contig_naming="ucsc")

    result = run_py(
        "import_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--id-vtable",
        id_vtable,
        "--output",
        out,
        "--genome-build",
        "GRCh37",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["chr1", "100", "rs1", "A", "G", ".", "0", "identity"]]
    meta_payload = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta_payload["target"]["genome_build"] == "GRCh38"
    assert meta_payload["target"]["contig_naming"] == "ucsc"


def test_import_sumstats_id_vtable_ignores_invalid_lookup_ids_and_audits_invalid_unmatched_and_ambiguous_source_ids(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    id_vtable = tmp_path / "lookup.vtable"
    out = tmp_path / "ss.vmap"
    write_lines(
        sumstats,
        [
            "SNP\tEA\tOA",
            ".\tA\tG",
            "rs_missing\tA\tG",
            "rs_dup\tA\tG",
            "rs_ok\tA\tG",
        ],
    )
    write_lines(meta, ["col_SNP: SNP", "col_EffectAllele: EA", "col_OtherAllele: OA"])
    write_vtable_with_meta(
        id_vtable,
        [
            "1\t10\t.\tA\tG",
            "1\t11\t\tA\tG",
            "1\t12\trs_dup\tA\tG",
            "1\t13\trs_dup\tA\tG",
            "1\t14\trs_ok\tA\tG",
        ],
        genome_build="GRCh37",
        contig_naming="ncbi",
    )

    result = run_py(
        "import_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--id-vtable",
        id_vtable,
        "--output",
        out,
    )

    assert result.returncode == 0, result.stderr
    assert "ignored 2 --id-vtable rows" in result.stderr
    assert read_tsv(out) == [["1", "14", "rs_ok", "A", "G", ".", "3", "identity"]]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "reason"],
        [".", "0", "invalid_id"],
        [".", "1", "id_not_found"],
        [".", "2", "ambiguous_id_match"],
    ]


def test_import_sumstats_id_vtable_requires_metadata_to_omit_chr_and_pos(tmp_path):
    sumstats = tmp_path / "ss.tsv"
    meta = tmp_path / "ss.yaml"
    id_vtable = tmp_path / "lookup.vtable"
    out = tmp_path / "ss.vmap"
    write_lines(sumstats, ["CHR\tPOS\tSNP\tEA\tOA", "1\t100\trs1\tA\tG"])
    write_lines(
        meta,
        [
            "col_CHR: CHR",
            "col_POS: POS",
            "col_SNP: SNP",
            "col_EffectAllele: EA",
            "col_OtherAllele: OA",
        ],
    )
    write_vtable_with_meta(id_vtable, ["1\t100\trs1\tA\tG"])

    result = run_py(
        "import_sumstats.py",
        "--input",
        sumstats,
        "--sumstats-metadata",
        meta,
        "--id-vtable",
        id_vtable,
        "--output",
        out,
    )

    assert result.returncode != 0
    assert "--id-vtable requires metadata to omit col_CHR and col_POS" in result.stderr

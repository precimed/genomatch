import json

from utils import read_tsv, run_py, write_json, write_lines


def test_normalize_contigs_vtable_ncbi_to_ucsc(tmp_path):
    source = tmp_path / "in.vtable"
    write_lines(source, ["1\t100\trs1\tA\tG", "MT\t200\trs2\tC\tT"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )
    out = tmp_path / "out.vtable"
    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "ucsc")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["chr1", "100", "rs1", "A", "G"], ["chrM", "200", "rs2", "C", "T"]]


def test_normalize_contigs_vtable_plink_to_ncbi(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["23\t100\trs1\tA\tG", "26\t200\trs2\tC\tT"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "plink"},
    )

    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "ncbi")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["X", "100", "rs1", "A", "G"], ["MT", "200", "rs2", "C", "T"]]


def test_normalize_contigs_vtable_plink_splitx_to_ncbi(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["25\t100\trs1\tA\tG", "24\t200\trs2\tC\tT"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "plink_splitx"},
    )

    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "ncbi")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["X", "100", "rs1", "A", "G"], ["Y", "200", "rs2", "C", "T"]]


def test_normalize_contigs_aliases_par1_and_par2_to_x(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["PAR1\t100\trs1\tA\tG", "PAR2\t200\trs2\tC\tT", "chrPAR1\t300\trs3\tG\tA"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37"},
    )

    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "ncbi")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["X", "100", "rs1", "A", "G"], ["X", "200", "rs2", "C", "T"], ["X", "300", "rs3", "G", "A"]]


def test_normalize_contigs_vtable_ncbi_to_plink_splitx_uses_grch37_par_boundaries(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(
        source,
        [
            "X\t60000\trs_before\tA\tG",
            "X\t60001\trs_par1_start\tA\tG",
            "X\t2699520\trs_par1_end\tA\tG",
            "X\t2699521\trs_after\tA\tG",
            "X\t154931043\trs_before_par2\tA\tG",
            "X\t154931044\trs_par2_start\tA\tG",
            "X\t155260560\trs_par2_end\tA\tG",
            "X\t155260561\trs_after_par2\tA\tG",
            "1\t100\trs_auto\tC\tT",
            "Y\t200\trs_y\tC\tT",
            "MT\t300\trs_mt\tC\tT",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "plink_splitx")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["23", "60000", "rs_before", "A", "G"],
        ["25", "60001", "rs_par1_start", "A", "G"],
        ["25", "2699520", "rs_par1_end", "A", "G"],
        ["23", "2699521", "rs_after", "A", "G"],
        ["23", "154931043", "rs_before_par2", "A", "G"],
        ["25", "154931044", "rs_par2_start", "A", "G"],
        ["25", "155260560", "rs_par2_end", "A", "G"],
        ["23", "155260561", "rs_after_par2", "A", "G"],
        ["1", "100", "rs_auto", "C", "T"],
        ["24", "200", "rs_y", "C", "T"],
        ["26", "300", "rs_mt", "C", "T"],
    ]
    out_meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert out_meta["contig_naming"] == "plink_splitx"
    assert out_meta["genome_build"] == "GRCh37"


def test_normalize_contigs_vtable_ncbi_to_plink_splitx_uses_grch38_par_boundaries(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(
        source,
        [
            "X\t10000\trs_before\tA\tG",
            "X\t10001\trs_par1_start\tA\tG",
            "X\t2781479\trs_par1_end\tA\tG",
            "X\t2781480\trs_after\tA\tG",
            "X\t155701382\trs_before_par2\tA\tG",
            "X\t155701383\trs_par2_start\tA\tG",
            "X\t156030895\trs_par2_end\tA\tG",
            "X\t156030896\trs_after_par2\tA\tG",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": "ncbi"},
    )

    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "plink_splitx")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["23", "10000", "rs_before", "A", "G"],
        ["25", "10001", "rs_par1_start", "A", "G"],
        ["25", "2781479", "rs_par1_end", "A", "G"],
        ["23", "2781480", "rs_after", "A", "G"],
        ["23", "155701382", "rs_before_par2", "A", "G"],
        ["25", "155701383", "rs_par2_start", "A", "G"],
        ["25", "156030895", "rs_par2_end", "A", "G"],
        ["23", "156030896", "rs_after_par2", "A", "G"],
    ]


def test_normalize_contigs_to_plink_splitx_requires_known_build(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["X\t60001\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "unknown", "contig_naming": "ncbi"},
    )

    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "plink_splitx")
    assert result.returncode != 0
    assert "requires metadata genome_build=GRCh37 or GRCh38" in result.stderr


def test_normalize_contigs_to_plink_splitx_rejects_unsupported_build(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["X\t60001\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "hg18", "contig_naming": "ncbi"},
    )

    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "plink_splitx")
    assert result.returncode != 0
    assert "unsupported genome_build in metadata" in result.stderr


def test_normalize_contigs_vmap_changes_target_side_only(tmp_path):
    source = tmp_path / "in.vmap"
    out = tmp_path / "out.vmap"
    write_lines(source, ["23\t100\tt1\tA\tG\t23\t3\tidentity", "26\t200\tt2\tC\tT\t.\t-1\tmissing"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {
            "object_type": "variant_map",
            "target": {"genome_build": "GRCh37", "contig_naming": "plink"},
        },
    )

    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "ncbi")
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["X", "100", "t1", "A", "G", "23", "3", "identity"], ["MT", "200", "t2", "C", "T", ".", "-1", "missing"]]


def test_normalize_contigs_vmap_drops_duplicate_targets_to_qc(tmp_path):
    source = tmp_path / "in.vmap"
    out = tmp_path / "out.vmap"
    write_lines(
        source,
        [
            "chr1\t100\tt1\tA\tG\tshardA\t0\tidentity",
            "1\t100\tt2\tA\tG\tshardB\t1\tidentity",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {
            "object_type": "variant_map",
            "target": {"genome_build": "GRCh37"},
        },
    )

    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "ncbi")
    assert result.returncode == 0, result.stderr
    assert "dropped 0 unresolved rows and 2 duplicate target rows" in result.stderr
    assert read_tsv(out) == []
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        ["shardA", "0", "t1", "duplicate_target"],
        ["shardB", "1", "t2", "duplicate_target"],
    ]


def test_normalize_contigs_drops_unresolved_rows_and_sets_output_naming(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["chr1\t100\trs1\tA\tG", "2\t200\trs2\tC\tT", "chrUn\t300\trs3\tG\tA"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37"},
    )

    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "ncbi")
    assert result.returncode == 0, result.stderr
    assert "dropped 1 unresolved rows" in result.stderr
    assert read_tsv(out) == [["1", "100", "rs1", "A", "G"], ["2", "200", "rs2", "C", "T"]]
    out_meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert out_meta["contig_naming"] == "ncbi"


def test_normalize_contigs_requires_metadata_sidecar(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["1\t100\trs1\tA\tG"])

    result = run_py("normalize_contigs.py", "--input", source, "--output", out, "--to", "ucsc")
    assert result.returncode != 0
    assert "metadata sidecar not found" in result.stderr.lower()

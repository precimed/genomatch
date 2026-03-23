import json

from utils import read_tsv, run_py, write_json, write_lines


def test_restrict_contigs_vtable_emits_vtable(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(
        source,
        [
            "1\t100\trs1\tA\tG",
            "2\t200\trs2\tC\tT",
            "2\t300\trs3\tG\tA",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("restrict_contigs.py", "--input", source, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "100", "rs1", "A", "G"],
        ["2", "200", "rs2", "C", "T"],
        ["2", "300", "rs3", "G", "A"],
    ]
    out_meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert out_meta == {
        "object_type": "variant_table",
        "genome_build": "GRCh37",
        "contig_naming": "ncbi",
    }


def test_restrict_contigs_vmap_preserves_original_provenance_and_type(tmp_path):
    source = tmp_path / "in.vmap"
    out = tmp_path / "out.vmap"
    write_lines(
        source,
        [
            "1\t100\tmid1\tA\tG\tshardA\t8\tswap",
            "1\t200\tmid2\tT\tC\tshardA\t3\tflip",
            "2\t300\tmid3\tC\tA\t.\t-1\tmissing",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {
            "object_type": "variant_map",
            "target": {"genome_build": "GRCh38", "contig_naming": "ncbi"},
        },
    )

    result = run_py("restrict_contigs.py", "--input", source, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "100", "mid1", "A", "G", "shardA", "8", "swap"],
        ["1", "200", "mid2", "T", "C", "shardA", "3", "flip"],
        ["2", "300", "mid3", "C", "A", ".", "-1", "missing"],
    ]
    out_meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert out_meta == {
        "object_type": "variant_map",
        "target": {"genome_build": "GRCh38", "contig_naming": "ncbi"},
    }


def test_restrict_contigs_filters_target_rows_by_chr2use(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(
        source,
        [
            "1\t100\trs1\tA\tG",
            "2\t200\trs2\tC\tT",
            "3\t300\trs3\tG\tA",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("restrict_contigs.py", "--input", source, "--output", out, "--contigs", "2")
    assert result.returncode == 0, result.stderr
    assert "outside --chr2use" in result.stderr
    assert read_tsv(out) == [["2", "200", "rs2", "C", "T"]]


def test_restrict_contigs_requires_declared_contig_naming(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["1\t100\trs1\tA\tG", "unknown\t200\trs2\tC\tT"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37"},
    )

    result = run_py("restrict_contigs.py", "--input", source, "--output", out)
    assert result.returncode != 0
    assert "contig_naming" in result.stderr
    assert "normalize_contigs.py" in result.stderr


def test_restrict_contigs_rejects_invalid_declared_contigs(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["chr1\t100\trs1\tA\tG", "unknown\t200\trs2\tC\tT"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("restrict_contigs.py", "--input", source, "--output", out)
    assert result.returncode != 0
    assert "inconsistent with declared contig_naming" in result.stderr


def test_restrict_contigs_rejects_unknown_contigs_under_declared_metadata(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["1\t100\trs1\tA\tG", "unknown\t200\trs2\tC\tT"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("restrict_contigs.py", "--input", source, "--output", out)
    assert result.returncode != 0
    assert "normalize_contigs.py" in result.stderr

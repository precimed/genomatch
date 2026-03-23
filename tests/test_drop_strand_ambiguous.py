import json

from utils import read_tsv, run_py, write_json, write_lines


def test_drop_strand_ambiguous_vtable_emits_vtable(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(
        source,
        [
            "1\t100\trs1\tA\tG",
            "1\t200\trs2\tA\tT",
            "2\t300\trs3\tAG\tCT",
            "2\t400\trs4\tA\tAT",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("drop_strand_ambiguous.py", "--input", source, "--output", out)
    assert result.returncode == 0, result.stderr
    assert "dropped 2 strand-ambiguous rows" in result.stderr
    assert read_tsv(out) == [
        ["1", "100", "rs1", "A", "G"],
        ["2", "400", "rs4", "A", "AT"],
    ]
    out_meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert out_meta == {
        "object_type": "variant_table",
        "genome_build": "GRCh37",
        "contig_naming": "ncbi",
    }


def test_drop_strand_ambiguous_vmap_preserves_original_provenance(tmp_path):
    source = tmp_path / "in.vmap"
    out = tmp_path / "out.vmap"
    write_lines(
        source,
        [
            "1\t100\tmid1\tA\tT\tshardA\t8\tidentity",
            "1\t200\tmid2\tA\tG\tshardA\t3\tswap",
            "2\t300\tmid3\tAG\tCT\tshardB\t5\tidentity",
        ],
    )
    write_json(
        source.with_name(source.name + ".meta.json"),
        {
            "object_type": "variant_map",
            "target": {"genome_build": "GRCh38", "contig_naming": "ncbi"},
        },
    )

    result = run_py("drop_strand_ambiguous.py", "--input", source, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "200", "mid2", "A", "G", "shardA", "3", "swap"]]


def test_drop_strand_ambiguous_fails_on_non_acgt_alleles(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["1\t100\trs1\tA\tN"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("drop_strand_ambiguous.py", "--input", source, "--output", out)
    assert result.returncode != 0
    assert "invalid allele code" in result.stderr.lower()


def test_drop_strand_ambiguous_requires_declared_contig_naming(tmp_path):
    source = tmp_path / "in.vtable"
    out = tmp_path / "out.vtable"
    write_lines(source, ["1\t100\trs1\tA\tG"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37"},
    )

    result = run_py("drop_strand_ambiguous.py", "--input", source, "--output", out)
    assert result.returncode != 0
    assert "contig_naming" in result.stderr

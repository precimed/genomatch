import pytest

from genomatch.vtable_utils import (
    VariantRow,
    read_vmap,
    read_vtable,
    validate_vmap_metadata,
    validate_vtable_metadata,
    write_vtable,
)


def test_read_vtable_rejects_invalid_pos(tmp_path):
    source = tmp_path / "bad.vtable"
    source.write_text("1\t0\trs1\tA\tG\n", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid pos"):
        read_vtable(source)


def test_write_vtable_rejects_invalid_rows(tmp_path):
    out = tmp_path / "bad.vtable"

    with pytest.raises(ValueError, match="missing id"):
        write_vtable(out, [VariantRow("1", "100", "", "A", "G")])


def test_read_vtable_rejects_non_acgt_alleles(tmp_path):
    source = tmp_path / "bad.vtable"
    source.write_text("1\t100\trs1\tA\tN\n", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid allele code"):
        read_vtable(source)


def test_read_vtable_normalizes_alleles_to_uppercase(tmp_path):
    source = tmp_path / "mixed_case.vtable"
    source.write_text("1\t100\trs1\ta\tc\n", encoding="utf-8")

    rows = read_vtable(source)
    assert rows == [VariantRow("1", "100", "rs1", "A", "C")]


def test_read_vmap_normalizes_alleles_to_uppercase(tmp_path):
    source = tmp_path / "mixed_case.vmap"
    source.write_text("1\t100\trs1\ta\tc\t.\t0\tidentity\n", encoding="utf-8")

    rows = read_vmap(source)
    assert rows[0].a1 == "A"
    assert rows[0].a2 == "C"


def test_metadata_variants_count_is_optional_but_validated_when_present():
    validate_vtable_metadata({"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    validate_vmap_metadata(
        {
            "object_type": "variant_map",
            "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"},
            "variants_count": 0,
        }
    )
    with pytest.raises(ValueError, match="variants_count"):
        validate_vtable_metadata(
            {
                "object_type": "variant_table",
                "genome_build": "GRCh37",
                "contig_naming": "ncbi",
                "variants_count": "1",
            }
        )

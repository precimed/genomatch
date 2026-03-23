import pytest

from genomatch.vtable_utils import VariantRow, read_vtable, write_vtable


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

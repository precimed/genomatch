import pandas as pd
import pytest

from genomatch.tabular_rows import VMapRowsTable, VariantRowsTable
from genomatch.vtable_utils import VMapRow, VariantRow


def test_variant_rows_table_roundtrip_rows() -> None:
    rows = [
        VariantRow("1", "101", "rs1", "A", "G"),
        VariantRow("2", "202", "rs2", "C", "T"),
    ]
    table = VariantRowsTable.from_rows(rows)
    assert list(table.frame.columns) == ["chrom", "pos", "id", "a1", "a2", "row_idx"]
    assert table.to_rows() == rows


def test_vmap_rows_table_roundtrip_rows() -> None:
    rows = [
        VMapRow("1", "101", "rs1", "A", "G", "shardA", 0, "identity"),
        VMapRow("2", "202", "rs2", "C", "T", "shardA", 1, "flip"),
    ]
    table = VMapRowsTable.from_rows(rows)
    assert list(table.frame.columns) == [
        "chrom",
        "pos",
        "id",
        "a1",
        "a2",
        "source_shard",
        "source_index",
        "allele_op",
        "row_idx",
    ]
    assert table.to_rows() == rows


def test_variant_rows_table_requires_columns() -> None:
    with pytest.raises(ValueError, match="variant frame is missing required columns: a2"):
        VariantRowsTable.from_frame(pd.DataFrame({"chrom": ["1"], "pos": ["1"], "id": ["rs1"], "a1": ["A"]}))


def test_vmap_rows_table_requires_columns() -> None:
    with pytest.raises(ValueError, match="vmap frame is missing required columns: allele_op"):
        VMapRowsTable.from_frame(
            pd.DataFrame(
                {
                    "chrom": ["1"],
                    "pos": ["1"],
                    "id": ["rs1"],
                    "a1": ["A"],
                    "a2": ["G"],
                    "source_shard": ["."],
                    "source_index": [0],
                }
            )
        )


def test_vmap_rows_table_to_rows_coerces_source_index_to_int() -> None:
    table = VMapRowsTable.from_frame(
        pd.DataFrame(
            {
                "chrom": ["1"],
                "pos": ["1"],
                "id": ["rs1"],
                "a1": ["A"],
                "a2": ["G"],
                "source_shard": ["."],
                "source_index": ["7"],
                "allele_op": ["swap"],
            }
        )
    )
    rows = table.to_rows()
    assert rows == [VMapRow("1", "1", "rs1", "A", "G", ".", 7, "swap")]

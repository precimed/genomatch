from genomatch.contig_cleanup_utils import normalize_target_rows
from genomatch.vtable_utils import VariantRow


def test_shared_contig_cleanup_helpers_align_unknown_marking_and_restriction():
    rows = [
        VariantRow("chr1", "100", "rs1", "A", "G"),
        VariantRow("chrUn", "200", "rs2", "C", "T"),
        VariantRow("26", "300", "rs3", "G", "A"),
    ]

    normalized = normalize_target_rows(rows, "ncbi")
    assert [row.chrom for row in normalized.rows] == ["1", "MT"]
    assert normalized.normalized_count == 2
    assert normalized.unknown_count == 1

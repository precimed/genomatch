import sys
from pathlib import Path

MATCH_DIR = Path(__file__).resolve().parents[1] / "match"
if str(MATCH_DIR) not in sys.path:
    sys.path.insert(0, str(MATCH_DIR))

from contig_cleanup_utils import normalize_target_rows
from vtable_utils import VariantRow


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

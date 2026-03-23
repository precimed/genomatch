import re
import sys
from pathlib import Path

import yaml

MATCH_DIR = Path(__file__).resolve().parents[1] / "match"
if str(MATCH_DIR) not in sys.path:
    sys.path.insert(0, str(MATCH_DIR))


def _schema_allows(schema: dict, value: object) -> bool:
    if "$ref" in schema:
        ref = schema["$ref"]
        assert isinstance(ref, str) and ref.startswith("#/definitions/")
        return _schema_allows(SCHEMA["definitions"][ref.rsplit("/", 1)[-1]], value)
    if "anyOf" in schema:
        return any(_schema_allows(option, value) for option in schema["anyOf"])
    schema_type = schema.get("type")
    if schema_type == "integer":
        if not isinstance(value, int) or isinstance(value, bool):
            return False
        minimum = schema.get("minimum")
        maximum = schema.get("maximum")
        return (minimum is None or value >= minimum) and (maximum is None or value <= maximum)
    if schema_type == "string":
        if not isinstance(value, str):
            return False
        pattern = schema.get("pattern")
        if pattern is None:
            return True
        return re.fullmatch(pattern, value) is not None
    raise AssertionError(f"unsupported schema fragment in test: {schema}")


SCHEMA = yaml.safe_load((MATCH_DIR / "schemas" / "cleaned-sumstats.yaml").read_text(encoding="utf-8"))
CHR_SCHEMA = SCHEMA["items"]["properties"]["CHR"]


def test_cleaned_sumstats_schema_accepts_supported_contig_namings():
    valid = [
        1,
        22,
        23,
        24,
        25,
        26,
        "1",
        "22",
        "23",
        "24",
        "25",
        "26",
        "X",
        "Y",
        "MT",
        "chr1",
        "chr22",
        "chrX",
        "chrY",
        "chrM",
    ]

    for value in valid:
        assert _schema_allows(CHR_SCHEMA, value), value


def test_cleaned_sumstats_schema_rejects_unsupported_contig_labels():
    invalid = [0, 27, "0", "27", "chr0", "chr23", "chrMT", "M", "XY", "chrXY", "unknown"]

    for value in invalid:
        assert not _schema_allows(CHR_SCHEMA, value), value

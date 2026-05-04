from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


SAMPLE_ID_MODE_FID_IID = "fid_iid"
SAMPLE_ID_MODE_IID = "iid"
SAMPLE_ID_MODE_CHOICES = (SAMPLE_ID_MODE_FID_IID, SAMPLE_ID_MODE_IID)


def normalize_sex_value(raw: str) -> int:
    try:
        sex = int(raw)
    except ValueError:
        return 0
    if sex not in (0, 1, 2):
        return 0
    return sex


def sex_to_label(sex: int) -> str:
    if sex == 1:
        return "1"
    if sex == 2:
        return "2"
    return "0"


def subject_key(*, fid: str, iid: str, sample_id_mode: str) -> Tuple[str, ...]:
    if sample_id_mode == SAMPLE_ID_MODE_IID:
        return (iid,)
    if sample_id_mode == SAMPLE_ID_MODE_FID_IID:
        return (fid, iid)
    raise ValueError(f"unsupported --sample-id-mode: {sample_id_mode!r}")


@dataclass(frozen=True)
class ParsedSampleRow:
    fid: str
    iid: str
    sex: int

    @property
    def display_id(self) -> str:
        return f"{self.fid}:{self.iid}"


@dataclass(frozen=True)
class ParsedSampleTable:
    path: Path
    kind: str
    rows: Tuple[ParsedSampleRow, ...]
    keys: Tuple[Tuple[str, ...], ...]
    key_to_index: Dict[Tuple[str, ...], int]
    signature: Tuple[object, ...]
    has_fid: bool | None = None

    @property
    def sample_count(self) -> int:
        return len(self.rows)

    @property
    def sexes(self) -> List[int]:
        return [row.sex for row in self.rows]

    @property
    def sample_ids(self) -> List[str]:
        return [row.display_id for row in self.rows]


@dataclass(frozen=True)
class SampleAxisPlan:
    output_sample_path: Path
    output_rows: Tuple[ParsedSampleRow, ...]
    source_tables: Dict[str, ParsedSampleTable]
    source_local_to_output: Dict[str, Tuple[int, ...]]
    source_present_output_indices: Dict[str, frozenset[int]]
    reconciliation_active: bool

    @property
    def output_sample_count(self) -> int:
        return len(self.output_rows)

    @property
    def output_sexes(self) -> List[int]:
        return [row.sex for row in self.output_rows]

    @property
    def output_sample_ids(self) -> List[str]:
        return [row.display_id for row in self.output_rows]


@dataclass(frozen=True)
class ReconciliationMissingnessSummary:
    mapped_variant_count: int
    output_sample_count: int
    total_missing_cells: int
    subjects_over_threshold: int
    variants_over_threshold: int


def _check_duplicate_key(
    key_to_index: Dict[Tuple[str, ...], int],
    key: Tuple[str, ...],
    *,
    path: Path,
    label: str,
) -> None:
    if key in key_to_index:
        raise ValueError(f"duplicate subject key in {label}: {path}")


def parse_fam_table(path: Path, *, sample_id_mode: str, label: str) -> ParsedSampleTable:
    rows: List[ParsedSampleRow] = []
    keys: List[Tuple[str, ...]] = []
    key_to_index: Dict[Tuple[str, ...], int] = {}
    signature: List[object] = []
    with open(path, "r", encoding="utf-8", newline="\n") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 6:
                raise ValueError(f"invalid .fam row in {path}")
            fid, iid, father, mother, sex_raw, pheno = parts[:6]
            row = ParsedSampleRow(fid=fid, iid=iid, sex=normalize_sex_value(sex_raw))
            key = subject_key(fid=fid, iid=iid, sample_id_mode=sample_id_mode)
            _check_duplicate_key(key_to_index, key, path=path, label=label)
            key_to_index[key] = len(rows)
            rows.append(row)
            keys.append(key)
            signature.append((fid, iid, father, mother, sex_raw, pheno))
    if not rows:
        raise ValueError(f"empty .fam: {path}")
    return ParsedSampleTable(
        path=path,
        kind="fam",
        rows=tuple(rows),
        keys=tuple(keys),
        key_to_index=key_to_index,
        signature=tuple(signature),
        has_fid=True,
    )


def parse_psam_table(path: Path, *, sample_id_mode: str, label: str) -> ParsedSampleTable:
    lines = [line.rstrip("\n") for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"empty .psam: {path}")
    header_parts = lines[0].split()
    if not header_parts or not header_parts[0].startswith("#"):
        raise ValueError(f"source .psam header not found: {path}")
    header_parts[0] = header_parts[0][1:]
    columns = {name: idx for idx, name in enumerate(header_parts)}
    if "IID" not in columns:
        raise ValueError(f".psam is missing required IID column: {path}")
    has_fid = "FID" in columns
    rows: List[ParsedSampleRow] = []
    keys: List[Tuple[str, ...]] = []
    key_to_index: Dict[Tuple[str, ...], int] = {}
    for line in lines[1:]:
        parts = line.split()
        max_idx = max(columns.values())
        if len(parts) <= max_idx:
            raise ValueError(f"malformed .psam row in {path}")
        fid = parts[columns["FID"]] if has_fid else "0"
        iid = parts[columns["IID"]]
        sex_raw = parts[columns["SEX"]] if "SEX" in columns else "0"
        row = ParsedSampleRow(fid=fid, iid=iid, sex=normalize_sex_value(sex_raw))
        key = subject_key(fid=fid, iid=iid, sample_id_mode=sample_id_mode)
        _check_duplicate_key(key_to_index, key, path=path, label=label)
        key_to_index[key] = len(rows)
        rows.append(row)
        keys.append(key)
    if not rows:
        raise ValueError(f"empty .psam sample body: {path}")
    return ParsedSampleTable(
        path=path,
        kind="psam",
        rows=tuple(rows),
        keys=tuple(keys),
        key_to_index=key_to_index,
        signature=tuple(lines),
        has_fid=has_fid,
    )


def require_psam_fid_presence_consistent(
    tables: Iterable[ParsedSampleTable],
    *,
    label: str,
) -> None:
    seen: set[bool] = set()
    for table in tables:
        if table.kind != "psam":
            continue
        assert table.has_fid is not None
        seen.add(bool(table.has_fid))
    if len(seen) > 1:
        raise ValueError(
            f"inconsistent .psam FID header presence under --sample-id-mode=fid_iid: {label}"
        )


def require_identical_sample_signatures(
    tables_by_shard: Dict[str, ParsedSampleTable],
    *,
    descriptor: str,
) -> Path:
    first_shard: str | None = None
    first_table: ParsedSampleTable | None = None
    for shard in sorted(tables_by_shard):
        table = tables_by_shard[shard]
        if first_table is None:
            first_shard = shard
            first_table = table
            continue
        if table.signature != first_table.signature:
            raise ValueError(
                f"source {descriptor} shard mismatch: all referenced source shards must have identical "
                f"{descriptor} contents; first mismatch at source_shard={shard!r}: {table.path}"
            )
    if first_table is None:
        raise ValueError(f"internal error: missing source {descriptor} table")
    assert first_shard is not None
    return first_table.path


def native_sample_axis_table_for_output_shard(
    tables_by_shard: Dict[str, ParsedSampleTable],
    source_shards: Sequence[str],
    *,
    descriptor: str,
    output_label: str,
) -> ParsedSampleTable:
    if not source_shards:
        raise ValueError(
            f"--sample-axis native cannot emit {output_label}: no mapped source rows define a native {descriptor} axis"
        )
    first_shard = source_shards[0]
    first_table = tables_by_shard[first_shard]
    for shard in source_shards[1:]:
        table = tables_by_shard[shard]
        if table.signature != first_table.signature:
            raise ValueError(
                f"--sample-axis native cannot emit {output_label}: mapped rows reference source shards "
                f"with different {descriptor} contents; first mismatch at source_shard={shard!r}: {table.path}"
            )
    return first_table


def mapped_source_shards_for_output_indices(indices: Sequence[int], vmap_rows: Sequence[object]) -> List[str]:
    seen: set[str] = set()
    out: List[str] = []
    for idx in indices:
        row = vmap_rows[idx]
        if int(row.source_index) == -1 or str(row.source_shard) in seen:
            continue
        shard = str(row.source_shard)
        seen.add(shard)
        out.append(shard)
    return out


def build_native_sample_axis_plan_for_output_shard(
    source_tables_by_shard: Dict[str, ParsedSampleTable],
    source_shards: Sequence[str],
    *,
    descriptor: str,
    output_label: str,
) -> SampleAxisPlan:
    native_table = native_sample_axis_table_for_output_shard(
        source_tables_by_shard,
        source_shards,
        descriptor=descriptor,
        output_label=output_label,
    )
    return build_sample_axis_plan(
        {shard: source_tables_by_shard[shard] for shard in source_shards},
        output_sample_path=native_table.path,
        explicit_target_table=None,
    )


def build_sample_axis_plan(
    source_tables_by_shard: Dict[str, ParsedSampleTable],
    *,
    output_sample_path: Path,
    explicit_target_table: ParsedSampleTable | None,
) -> SampleAxisPlan:
    if explicit_target_table is None:
        output_table = next(iter(source_tables_by_shard.values()))
        source_local_to_output = {
            shard: tuple(range(table.sample_count))
            for shard, table in source_tables_by_shard.items()
        }
        source_present_output_indices = {
            shard: frozenset(range(table.sample_count))
            for shard, table in source_tables_by_shard.items()
        }
        return SampleAxisPlan(
            output_sample_path=output_sample_path,
            output_rows=output_table.rows,
            source_tables=source_tables_by_shard,
            source_local_to_output=source_local_to_output,
            source_present_output_indices=source_present_output_indices,
            reconciliation_active=False,
        )

    source_local_to_output: Dict[str, Tuple[int, ...]] = {}
    source_present_output_indices: Dict[str, frozenset[int]] = {}
    for shard, table in source_tables_by_shard.items():
        local_to_output = tuple(
            explicit_target_table.key_to_index.get(key, -1)
            for key in table.keys
        )
        source_local_to_output[shard] = local_to_output
        source_present_output_indices[shard] = frozenset(idx for idx in local_to_output if idx != -1)
    return SampleAxisPlan(
        output_sample_path=output_sample_path,
        output_rows=explicit_target_table.rows,
        source_tables=source_tables_by_shard,
        source_local_to_output=source_local_to_output,
        source_present_output_indices=source_present_output_indices,
        reconciliation_active=True,
    )


def compute_reconciliation_missingness_summary(
    plan: SampleAxisPlan,
    vmap_rows: Sequence[object],
) -> ReconciliationMissingnessSummary | None:
    if not plan.reconciliation_active:
        return None
    mapped_rows_by_shard: Dict[str, int] = {}
    for row in vmap_rows:
        if int(row.source_index) == -1:
            continue
        mapped_rows_by_shard[str(row.source_shard)] = mapped_rows_by_shard.get(str(row.source_shard), 0) + 1
    mapped_variant_count = sum(mapped_rows_by_shard.values())
    output_sample_count = plan.output_sample_count
    subject_missing_counts = [0] * output_sample_count
    total_missing_cells = 0
    variants_over_threshold = 0
    for shard, mapped_count in mapped_rows_by_shard.items():
        present = plan.source_present_output_indices[shard]
        absent_count = output_sample_count - len(present)
        total_missing_cells += mapped_count * absent_count
        if absent_count * 2 > output_sample_count:
            variants_over_threshold += mapped_count
        if mapped_count == 0 or absent_count == 0:
            continue
        for output_idx in range(output_sample_count):
            if output_idx not in present:
                subject_missing_counts[output_idx] += mapped_count
    subjects_over_threshold = 0
    if mapped_variant_count:
        subjects_over_threshold = sum(
            1
            for missing_count in subject_missing_counts
            if missing_count * 2 > mapped_variant_count
        )
    return ReconciliationMissingnessSummary(
        mapped_variant_count=mapped_variant_count,
        output_sample_count=output_sample_count,
        total_missing_cells=total_missing_cells,
        subjects_over_threshold=subjects_over_threshold,
        variants_over_threshold=variants_over_threshold,
    )

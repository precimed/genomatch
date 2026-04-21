# Telemetry

## Scope

Telemetry is operational user feedback emitted during CLI execution.

This spec defines only high-level behavior. Message templates, exact wording, and logging style are implementation-defined.

Detailed implementation guidance lives in [IMPLEMENTATION_GUIDE.md](../IMPLEMENTATION_GUIDE.md).

## Contract

- Telemetry is emitted to `stderr`.
- Telemetry is non-API output and may change between versions.
- Telemetry must use Python `logging`.
- Logging configuration should include timestamps for every emitted message.
- Direct `print(..., file=sys.stderr)` is not allowed for telemetry.
- Wrappers (`prepare_variants.py`, `project_payload.py`) log invoked subcommands and wrapper status to `stderr`.
- Wrappers may forward child-process stdout/stderr to `stderr` during orchestration.
- Error handling remains defined by CLI behavior contracts (`exit code != 0` and user-visible failure reason).
- Caller/tool entrypoints own INFO-level start/progress/completion messaging.
- Utility functions should not emit INFO for normal flow; they should return counts/status for callers to report.
- Utility functions may emit WARNING/ERROR only for local anomalies that cannot be better contextualized by the caller.
- Avoid duplicate messaging for the same event across utility and caller boundaries.

## Levels and granularity

- Implementations should include enough context for users to understand progress, warnings, and failures.
- Typical categories are informational progress, recoverable warnings, and fatal errors.
- Per-tool telemetry shape is intentionally not standardized in this spec.

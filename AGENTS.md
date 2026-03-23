# AGENTS.md

## Source of truth

Follow SPEC.md and all documents in spec/ folder, README.md (user-facing documentation), TESTS.md, INSTALL.md as the source of truth. Do not invent undocumented behavior. When docs and code disagree, align implementation to the spec unless explicitly instructed otherwise.
When making changes to source-of-truth files, avoid repeating the information in multiple documents.

Within `README.md`, the following sections are intentionally normative and may serve as source of truth:

- the top-level object and metadata summary
- `Canonical tools`
- `Using .vmap vs .vtable as input across tools`

## Change discipline

Prefer minimal, surgical edits. Avoid unrelated refactors, formatting churn, and public-interface changes unless required.
Prefer small pure functions for transformation logic. Keep I/O separate from core logic. Preserve existing module organization.

## Testing discipline

Run the narrowest relevant tests first. Add or update tests for behavior changes. Never claim tests passed unless they were actually run. If tests could not be run, say so clearly.
Use `match-liftover` conda environment to run tests (non-trivial dependencies)

## CLI, schema stability, documentation

Preserve stable CLI flags, file formats, column semantics, schema fields unless explicitly instructed otherwise.
If CLI behavior, install steps, reference-data requirements, schemas, or examples change, update the corresponding docs in the same change.

## Domain guardrails

Fail early on ambiguous genome build, contig naming, allele conventions, or metadata. Do not silently drop unmatched variants; report counts and reasons. Keep liftover, contig normalization, and allele reconciliation auditable.

## Reproducibility and environment hygiene

Do not hardcode local paths, usernames, or machine-specific assumptions. Do not assume internet access at runtime. Keep new scripts runnable from a clean checkout after following INSTALL.md.

## Review checklist

Check edge cases around missing metadata, duplicate IDs, malformed input, backward compatibility, docs/test sync, and portability before finalizing.

## Github commits

Prefer informative commit message and detailed description.
Commit one logical unit of changes at a time.
For code changes, run the full test suite before commits. Not applicable to spec/docs change.

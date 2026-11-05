# Schema compatibility fixtures

These committed files were written by the initial schema `1.0` implementation.
They cover stable row order and IDs, multibase alleles, and, for `.vmap`, named
source shards, nonzero source indices, and non-identity allele operations.

Writer provenance: checkout based on commit `04f74dfcdb2bd2f478ae09d3d67c9e27b7f01f20`, with the initial
schema-version implementation applied (unreleased when generated).

Reproduce from the repository root in `match-liftover` using that writer checkout:

```sh
python tests/fixtures/schema/generate_1_0.py
```

The generator uses normal object and metadata writers and checks that their
current schema is `1.0`. It writes both objects and their `.meta.json` sidecars.
It is not run by tests or routine fixture regeneration. Preserve these bytes
when introducing later schemas; use the original writer to reproduce them.
Expected logical contents are asserted independently in `tests/test_schema_metadata.py`.

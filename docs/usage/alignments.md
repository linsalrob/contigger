# Alignment and evidence commands

`contigger merge` invokes minimap2 for planned candidate pairs. It batches queries by one target and validates/reuses target indexes. `--minimap2-preset` selects the assembly preset; `--index-dir` controls the cache.

The separate command:

```bash
contigger benchmark-junction-remapping --dataset test_data --preset map-ont
```

is a checked-in evidence benchmark. It remaps dataset-specific targeted reads and scores spanning observations. It is not a general raw-read merge interface, and its observations do not automatically authorize a production merge.

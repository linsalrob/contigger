# Contigger compact validation dataset

Version 1.0.0; seed 47291. This benchmark tests conservative exact, containment, overlap, threshold, circular, mutation, negative and repeat-ambiguity decisions under the rule **a missed merge is preferable to a false merge**.

It contains only benchmark-derived reference slices and targeted ONT reads. It excludes the complete reads, assembly, GFA, source BAM and indexes. Coordinates are zero-based, half-open; strand and circular wrapping are explicit in `expected/`. Use `manifest.tsv` with reads/BAMs or `manifest_no_reads.tsv` for contig-only runs. PAFs are realistic unfiltered asm5/asm20 all-vs-all classifier inputs. Source-contig BAMs support pileup and extraction; a newly merged sequence must be remapped before direct junction validation.

After unpacking beside `benchmarking/`, run `python benchmarking/validate_benchmark.py test_data`.

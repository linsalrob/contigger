# Benchmarks

Checked-in Pseudomonas-derived PAF and truth files test catalogue identity, candidate recall, relationship classification, graph ambiguity, policy, paths, and junction evidence. Baselines live under `benchmarks/` and are validated by CI; they are not claims of production sensitivity.

The merge regression preserves zero false constructed joins, including both directions of the known `end_tolerance_51` classifier-like case. The end-to-end sequence baseline is conservative: exact/RC aliases and eligible containments are represented, while unsupported overlap paths remain deferred.

`benchmarking/merge_scale.py` provides deterministic pre-alignment scale measurements. Current smoke measurements cover 10,000 contigs and a 100,000-contig short-sequence case. Record hardware, sequence lengths, candidate counts, peak RSS, and settings with any new measurement. Never update a baseline silently; explain every changed count and add a regression for every new false join.

For an external experimental collection, use `benchmarking/real_fasta_scale.py` to make
nested, seeded identifier-hash subsets and `benchmarking/profile_fasta_candidates.py` to
exercise the real catalogue/candidate path.  The helper reports its measured source
record and nucleotide-base totals; these should be preferred to filesystem-size-derived
estimates.  Generated FASTA files, manifests containing private paths, raw scheduler
logs, and full profiling JSON stay outside version control.  Commit only aggregate,
non-sensitive findings with their command settings and hardware context.

The checked-in aggregate real-contig candidate baseline is
`benchmarks/scale-results/real-contig-candidate-scaling.json`.
It records 10k, 100k, and 1m nested samples, but is explicitly pre-alignment and not a
claim that an end-to-end full-size merge will fit the same resources.

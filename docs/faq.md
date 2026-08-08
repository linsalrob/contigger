# FAQ

**Does Contigger replace an assembler?** No. It reconciles assembled contigs.

**Does it polish contigs?** No. It does not perform general polishing or variant calling.

**Why did two 99% identical contigs not merge?** Identity is only one condition; geometry, branches, repeats, conflicts, and policy also matter.

**Can I use only FASTA files?** Yes. This is the default sequence-only workflow.

**Do I need raw reads?** No for sequence-only merging. Raw FASTQ is not a direct `merge` input.

**Can I use ONT, Illumina, or HiFi reads?** Their mapped BAM/CRAM files can be validated if references match. No current reviewed merge policy resolves imperfect consensus for any technology automatically.

**Can I combine technologies?** You may label samples, but evidence remains sample-scoped and is not blindly pooled.

**Does merge use GFA graphs?** The manifest accepts and validates `assembly_graph`, but graph edges do not currently authorize merges.

**Does it support circular contigs?** Treat circular-origin relationships as requiring review; no general circular consensus should be assumed.

**Why are output IDs different?** Output IDs are deterministic representatives or path IDs; provenance maps them back to every source identifier.

**Should I delete original assemblies?** No. Keep them with the manifest and evidence files.

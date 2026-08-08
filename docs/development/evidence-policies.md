# Evidence policies

`BamEvidenceProvider` validates per-sample alignment indexes, integrity, exact reference names/lengths, and captures samtools commands. Pileups remain source-contig evidence. Targeted remapping uses a provisional junction and bounded source-end reads in the evidence benchmark modules.

The policy model separates junction support from base/indel choice. Review artifacts carry truth and candidate-baseline digests, technology/preset, thresholds, and reviewer state. Current checked-in policies are not approved for production consensus. Consequently `--evidence alignments` validates and reports, but does not authorize imperfect overlap construction. This is intentional: a spanning observation or pooled majority cannot hide sample-specific variation.

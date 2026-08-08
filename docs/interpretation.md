# Interpreting results

A merged sequence means that the accepted path had compatible terminal geometry and an explicitly reconcilable overlap under the selected evidence mode. It does not prove that the source contigs came from the same biological molecule or strain.

Exact/RC deduplication is sequence identity, not biological identity. Containment is safer than approximate overlap because the surviving sequence is already present in full. Deferred joins may reflect repeats, strain variation, conflicting alleles, cycles, circular-origin ambiguity, or insufficient evidence. Conserved genes can create internal similarity without a valid join; mosaic phage genomes are especially likely to need manual review.

Retain original FASTA, BAM/CRAM, manifest, command line, and output stats. Inspect provenance for source coordinates and ambiguity for unresolved alternatives before downstream analysis.

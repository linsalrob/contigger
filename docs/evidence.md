# Read evidence

Contigger keeps four questions separate:

```text
sequence overlap → source BAM/CRAM evidence → targeted junction evidence → sequence decision
```

Source BAMs can provide sample-scoped depth, allele counts, qualities, and mapping context on original contigs. They cannot directly prove a junction that did not exist in the source reference. Targeted remapping to a provisional junction is a separate evidence operation.

The current merge supports `--evidence none` and `--evidence alignments`. The latter validates each sample's BAM/CRAM and writes `join_support.tsv`; imperfect overlaps remain `DEFERRED` because no checked-in reviewed policy authorizes choosing a SNP or indel. Contradictory sample-specific alleles are never pooled. `--evidence reads` is retained as a legacy CLI value but rejected by merge.

No read mode should be described as polishing, variant calling, or proof of biological identity. See [the developer evidence-policy notes](development/evidence-policies.md) for review artifacts and benchmark boundaries.

# Concepts

**Exact duplicate:** two source sequences have the same bases. **Reverse-complement duplicate:** one is the reverse complement of the other. The catalogue emits one representative while provenance retains both sources.

**Containment:** one contig lies end-to-end inside another. A uniquely eligible contained contig may be omitted from representative FASTA, but remains in provenance.

**Terminal overlap:** the suffix of one oriented contig matches the prefix of another. Internal similarity is not a join.

```text
Exact:       A ====================
             B ====================
Containment: A ==========================
             B      ============
Terminal:    A ==================||||
             B              ||||=================
Internal:    A ========||||||||========
             B     ====||||||||====
```

Repeats, strain variation, assembly errors, and circular-origin geometry can create plausible but ambiguous relationships. Contigger preserves those cases. Identity is a filter, not a consensus rule: 98% identity alone never authorizes a merge.

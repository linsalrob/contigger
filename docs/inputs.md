# Input files

## FASTA

Each sample must provide a readable FASTA. Contig identifiers must be unique within that file; sequence lines may contain whitespace and standard IUPAC DNA symbols (`ACGTRYSWKMBDHVN`). Empty sequences, duplicate identifiers, malformed headers, and invalid symbols fail validation. Plain and gzip-compressed text are supported by the FASTA reader; use `.fa`, `.fna`, or `.fasta` as appropriate.

## BAM/CRAM

An optional BAM or CRAM must have an adjacent BAI index for BAM or CRAI index for CRAM. For `validate-alignments` and `--evidence alignments`, it must be coordinate-readable and its reference names and lengths must exactly match the corresponding source FASTA. A source BAM supports sample-scoped observations on source contigs; it cannot, by itself, prove a newly constructed junction.

## Assembly graphs and raw reads

`assembly_graph` is accepted and file-checked in the manifest, but the merge pipeline does not currently use GFA graph edges to authorize merges. Raw FASTQ files are not manifest inputs to `contigger merge`. The checked-in junction-remapping benchmark command operates on its dataset rather than arbitrary raw-read paths.

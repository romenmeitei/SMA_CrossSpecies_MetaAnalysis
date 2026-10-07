# Publication gates

The analysis is not submission-ready until every required gate is satisfied.

## Required before any primary results are interpreted

- [ ] Manifest validation has zero errors.
- [ ] Sample annotations were checked against each GEO record and primary paper.
- [ ] Raw/processed input checksums and retrieval dates are recorded.
- [ ] Every gene-level table uses the same disease-minus-control orientation.
- [ ] Ensembl release and one-to-one ortholog file are archived.
- [ ] QC plots show no unresolved swaps, severe outliers, or design confounding.
- [ ] GSE207890 primary and mouse-level models both ran successfully.
- [ ] Single-cell evidence is described as localization, not replicated DE.
- [ ] Tests pass in a clean environment.

## Required for a strong journal submission

- [ ] RNA-seq datasets are reprocessed from FASTQ with the same genome release and
      quantification workflow, or processed-summary dependence is justified and
      tested as a sensitivity analysis.
- [ ] At least one additional independent human cohort/model is identified or the
      limited human replication is explicit in title, abstract, and limitations.
- [ ] Leave-one-dataset-out, strict-vs-sensitivity, heterogeneity, and tissue-stratum
      results are reported.
- [ ] Enrichment uses an explicit measured-gene universe and versioned libraries.
- [ ] Network claims survive degree-matched permutation on the full versioned graph.
- [ ] Candidate genes/pathways have orthogonal experimental or clinical validation.
- [ ] Code, parameter files, full results, and reporting checklist are archived with
      a permanent DOI.

## Claims that are prohibited by this design alone

- “Conserved” based only on symbol overlap.
- “Treatment reversed” without the treated-minus-disease contrast.
- “Replicated single-cell DE” from one pooled sample per condition.
- “Hub gene” from a graph assembled in disconnected API chunks.
- Causal, diagnostic, or therapeutic-efficacy claims without external validation.


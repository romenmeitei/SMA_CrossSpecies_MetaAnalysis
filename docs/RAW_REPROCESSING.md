# Confirmatory raw-read reprocessing profile

The Colab-feasible profile uses compact GEO supplementary files for early iteration.
Before submission, reprocess GSE108094, GSE69175, GSE56284, and GSE209926 from
FASTQ in one workflow.

## Locked recommendations

1. Retrieve run accessions from the NCBI SRA Run Selector and archive the run table.
2. Use one containerized workflow release (for example nf-core/rnaseq) with explicit
   genome assets for GRCh38 and GRCm39 from the same Ensembl release used for
   orthology.
3. Apply the same adapter/quality filters and transcript-to-gene aggregation rules.
4. Preserve donor/cell-line/animal design variables. Do not merge technical
   replicates until their identity has been confirmed.
5. Quantify genes with Ensembl stable IDs and preserve versioned IDs in the raw
   output; strip versions only in the integration table.
6. Run MultiQC and archive the report. Define exclusion rules before inspecting
   disease effects.
7. Fit counts with DESeq2 or edgeR and export the standardized result schema used by
   `results/de/`.

Raw-read reprocessing is deliberately not launched by the default Colab notebook:
it is compute- and storage-intensive and should run on an HPC/cloud workflow engine.
The downstream integration remains identical once standardized DE files are placed
in `results/de/`.


# SMA cross-species, tissue-aware and treatment-reversal meta-analysis

## Manuscript-aligned v1.0.0 release

This repository accompanies the cross-species, tissue-aware SMA
transcriptomic study and its frozen reproducibility release.

The primary small-k Hartung-Knapp and modified Hartung-Knapp analyses
produced no transcriptome-wide FDR-significant genes. Corrected REML
results are therefore interpreted as screening evidence.

FAM107A is the sole Tier-4 internally robust candidate. Its external
support is directionally concordant but not statistically confirmatory.
UGP2, Brain Development splicing enrichment and AMPA-receptor
trafficking are reported as secondary findings.

Complete outputs are distributed as checksummed GitHub Release assets.
See `docs/GITHUB_REPRODUCIBILITY.md` and `docs/RELEASE_ASSETS.md`.


This repository is a from-scratch replacement for the earlier “conserved-gene”
and chunked-network workflow. It preserves the biological question—what molecular
programs recur across SMA models and are reversed by treatment—but changes the
unit of evidence from gene-list overlap to study-level effect estimates.

No result is embedded or preordained. The repository downloads public data,
records checksums, validates study eligibility, estimates within-dataset effects,
maps only versioned one-to-one orthologs, performs random-effects meta-analysis,
and keeps neural discovery, skeletal-muscle reversal, and single-cell localization
as separate evidence layers.

## Start in Google Colab

1. Download the complete project ZIP.
2. Open `notebooks/00_Master_Run.ipynb` in Colab.
3. Run the first configuration cell. If the project is not already present, it
   will prompt you to upload the complete ZIP and will extract and locate it.
4. Continue from top to bottom. Heavy and optional stages are controlled
   by explicit switches; the default first run performs validation and small
   processed-data downloads only.

The master notebook is the auditable runner. The numbered notebooks are focused
companions for QC, DE, meta-analysis, treatment reversal, splicing, enrichment,
and final validation.

## Evidence architecture

| Layer | Datasets | Primary purpose | Key restriction |
|---|---|---|---|
| Human motor neuron | GSE108094 | Human discovery and splicing | Cell line is the unit; technical replicates are not independent |
| Human motor neuron | GSE69175 | Directional sensitivity | One control and one SMA line; disease and line are confounded |
| Mouse motor neuron | GSE56284 | Cross-species neural discovery | In-vitro mESC-derived motor neurons |
| Mouse spinal cord | GSE209926 | Intact-tissue neural validation | Bulk signal can include cell-composition change |
| Mouse myofiber | GSE207890 | Treatment reversal | Fibers are repeated measures within seven mice |
| Mouse single cell | GSE208629 | Cell-type localization only | One pooled sample per condition; cells are not replicates |

The official GEO records are linked in `config/dataset_manifest.csv`. GSE38417 is
explicitly prohibited because it is a Duchenne muscular dystrophy dataset, not an
SMA dataset. GSE154106 is excluded from transcriptomic DE because it is a
ribosome-profiling assay.

## Statistical design

### Within-dataset effects

- RNA-seq count matrices use DESeq2, with disease minus control as the invariant
  effect direction.
- Processed expression matrices use limma, with donor/cell-line blocking when
  repeated technical or differentiation replicates exist.
- GSE207890 uses `duplicateCorrelation` with mouse as the block. A prespecified
  sensitivity analysis averages fibers to one profile per mouse before fitting.
- Published Cuffdiff summaries can seed the Colab-feasible analysis. For a final
  submission, the raw-read harmonized profile described in
  `docs/RAW_REPROCESSING.md` is the preferred confirmatory analysis.

Every standardized result contains accession, gene ID, effect, standard error
when available, signed z statistic, p value, FDR, species, tissue, role,
eligibility, and contrast.

### Cross-species meta-analysis

- Human and mouse genes are linked only by a pinned Ensembl Compara
  `ortholog_one2one` table. Gene-symbol equality is never used as evidence of
  orthology.
- Motor-neuron effects are the primary cross-species stratum. Spinal cord is
  reported separately; a joint “neural” analysis is explicitly labeled a
  sensitivity model.
- When study standard errors are available, per-gene effects use random-effects
  REML and report pooled effect, 95% CI, tau-squared, Cochran Q, Q p value, and
  I-squared. If a source exposes only signed significance statistics, weighted
  Stouffer meta-analysis is used and no pooled fold change is claimed.
- Strict analyses exclude `sensitivity_only` studies. Leave-one-dataset-out and
  species concordance outputs are generated separately.

### Treatment reversal

GSE207890 supplies two independent contrasts:

1. `disease - control`
2. `treated - disease`

The continuous reversal score is `-z_disease × z_treatment`. A gene is called
restored only when directions oppose, both contrasts pass their prespecified FDR,
both absolute z scores pass the minimum, and the mouse-level sensitivity analysis
supports the result. This prevents “reversal” from being inferred from the disease
contrast alone.

### Splicing, enrichment, and networks

- Splicing is event-level: rMATS FDR and absolute delta-PSI are thresholded before
  gene aggregation. Gene names containing a substring are never treated as a
  biological module.
- Enrichment is preranked and records the tested gene-set library and ranking.
- Network analysis is optional and off by default. It consumes a complete,
  versioned STRING bulk graph and uses degree-matched permutations. API batching
  is prohibited because it can remove edges between batches.

The full prespecified analysis is in `docs/STATISTICAL_ANALYSIS_PLAN.md`.

## Repository map

```text
config/       locked dataset manifest, thresholds, and sample annotations
notebooks/    Colab master runner plus seven focused stage notebooks
src/          tested Python package for validation and evidence integration
R/            Bioconductor engines for DE and orthology export
scripts/      input preparation and notebook generation helpers
tests/        synthetic unit tests (never scientific results)
data/         downloaded, intermediate, processed, and audit files
results/      machine-readable analysis outputs
figures/      publication figures
docs/         SAP, raw-read profile, data dictionary, and publication gates
```

## Command-line execution

```bash
python -m pip install -e '.[enrichment,notebook,test]'
sma-xspecies validate
sma-xspecies download
sma-xspecies standardize
Rscript environment/install_bioc.R
Rscript R/fetch_orthologs.R resources/ensembl_human_mouse_one_to_one.tsv 113
pytest
```

Dataset-specific R calls and later stages are shown in the notebooks. Each run
writes validation/audit tables and `results/run_record.json`; scientific outputs
are not silently reused.

## Interpretation boundaries

This computational project can support prioritization and mechanistic hypotheses;
it cannot by itself establish therapeutic efficacy or causality. A publishable
manuscript should distinguish discovery, replication, localization, and treatment
reversal, report all sensitivity analyses, and include independent biological or
clinical validation for strong therapeutic claims.

## Primary data and methodological context

- [GSE108094: human SMA and control motor neurons](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108094)
- [GSE69175: purified patient-derived motor neurons](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69175)
- [GSE56284: mouse ESC-derived motor neurons](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56284)
- [GSE209926: severe-SMA mouse spinal cord](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209926)
- [GSE207890: SMA mouse single-myofiber and AFS treatment study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207890)
- [GSE208629: severe-SMA mouse spinal-cord single-cell atlas](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208629)
- [Prior comparative SMA transcriptomic meta-analysis](https://doi.org/10.1186/s12920-024-02040-0)
- [STRING database 2025 update](https://doi.org/10.1093/nar/gkae1113)

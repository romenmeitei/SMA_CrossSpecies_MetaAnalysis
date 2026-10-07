# Prespecified statistical analysis plan

Version: 0.1.0  
Locked parameters: `config/analysis_parameters.yaml`

## Objectives and estimands

### Primary objective

Estimate the direction and magnitude of SMA-associated expression change among
one-to-one human–mouse orthologs in motor-neuron models. The primary estimand is a
random-effects mean log2 fold change for `disease - control`. A gene requires at
least two eligible studies and representation from both species.

### Secondary objectives

1. Determine whether primary motor-neuron findings replicate in mouse spinal cord.
2. Quantify treatment reversal in mouse skeletal myofibers using two independently
   fitted contrasts.
3. Localize prioritized programs to cell types using GSE208629 without using cells
   as biological replicates.
4. Test event-level differential splicing in GSE108094.
5. Evaluate pathways by preranked GSEA and network connectivity by degree-matched
   permutation only after gene-level inference is complete.

## Eligibility and exclusions

Eligibility is encoded before analysis. Core datasets must have at least two
independent biological units per group, an admissible transcriptomic assay, a clear
SMA/control contrast, and sample annotations sufficient to identify repeated
measures. `sensitivity_only` and `localization_only` datasets cannot enter strict
primary inference. Assay, tissue, or disease mismatch triggers exclusion.

## Preprocessing and quality control

- Record source URL, retrieval status, byte size, and SHA-256 checksum.
- Verify sample IDs against the locked sample annotation.
- Inspect library size/distribution, missingness, PCA/MDS, sample correlation, and
  condition/batch confounding before fitting.
- RNA-seq counts: filter low expression and fit DESeq2. No quantile normalization.
- Microarray: log2 transform only if distributional diagnostics indicate raw-scale
  intensities; quantile-normalize, annotate, then choose one probe per gene by
  highest IQR before modeling.
- The analytical unit must be donor, cell line, animal, or independent culture as
  declared in the manifest. Technical replicates are blocked or aggregated.

## Per-study models

All effects are oriented as disease minus control. Treatment is treated minus
disease. The primary microarray model for GSE207890 is

`expression ~ 0 + condition`, with mouse as a duplicate-correlation block.

The sensitivity model averages expression across fibers within mouse and fits mice
as samples. Results unsupported by the animal-level model are not called robust.

## Orthology

The primary mapping is Ensembl release 113, human–mouse `ortholog_one2one` only.
Version suffixes are removed from Ensembl IDs. One-to-many and many-to-many
relationships are excluded from primary inference and may be explored only in a
separately labeled sensitivity analysis.

## Meta-analysis

For gene g and study i, the inverse-variance model is

`y_gi ~ Normal(mu_g, se_gi^2 + tau_g^2)`.

Tau-squared is estimated by REML. Report k, pooled effect, standard error, 95% CI,
z, p, FDR, tau-squared, Q, Q p value, I-squared, participating studies/species, and
sign concordance. When standard errors are unavailable, combine signed z scores by
weighted Stouffer and do not report a pooled fold change.

Benjamini–Hochberg FDR is applied within each declared analysis family: primary
motor-neuron meta-analysis, spinal-cord validation, joint-neural sensitivity,
treatment reversal, splicing event type, and each enrichment collection.

## Heterogeneity and sensitivity analyses

- Flag I-squared above 75%.
- Leave one dataset out when at least three studies contribute.
- Compare strict results with inclusion of GSE69175.
- Report motor neuron and spinal cord separately; joint neural results are secondary.
- Compare blocked-fiber and mouse-level GSE207890 estimates.
- Report human–mouse signed-z Spearman correlation and sign agreement.
- Evaluate whether results depend on minimum effect, FDR, and orthology restrictions.

## Treatment reversal

For each gene, calculate `R = -z_disease × z_treatment`. Positive values favor
opposing effects. A restored call requires: opposing effects; FDR at or below 0.10
for both contrasts; absolute z at least 1 for both; and restoration fraction between
0 and 1.5. Values above 1.5 are labeled overshoot. The same rule is applied to the
mouse-level sensitivity model; robust calls require both models.

## Reporting

Report exact sample counts and biological units, exclusions, missing inputs,
software versions, parameter file, input checksums, all tested genes/events, and
negative/null findings. Do not describe associations as causal or therapeutic
validation without independent experiments.


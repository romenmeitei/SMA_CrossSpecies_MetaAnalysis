# Output data dictionary

## Standardized differential-expression files

| Column | Definition |
|---|---|
| accession | GEO study accession |
| gene_id | Version-stripped Ensembl gene ID where available |
| gene_symbol | Display symbol; never used alone for orthology |
| effect | Disease-minus-control or treated-minus-disease log2 fold change |
| se | Standard error of the effect, if available |
| z | Signed study statistic; derived from p and direction only when SE/stat is absent |
| p_value | Nominal two-sided p value |
| fdr | Within-study Benjamini–Hochberg or source q value |
| species | Homo_sapiens or Mus_musculus |
| tissue | Prespecified tissue stratum |
| eligibility | eligible, sensitivity_only, or localization_only |
| contrast | Invariant contrast label |

## Gene meta-analysis

`estimate`, `se`, and confidence limits are emitted only for inverse-variance REML.
`tau2`, `q`, `q_p_value`, and `i2` quantify between-study heterogeneity. `datasets`
and `species` make the evidence set inspectable for every gene.

## Treatment reversal

`reversal_score` is the negative product of the disease and treatment z statistics.
`restoration_fraction` is `-treatment_effect / disease_effect`. `passes_conjunction`
requires both contrasts to meet their thresholds. `robust_reversal` additionally
requires the mouse-level sensitivity model.


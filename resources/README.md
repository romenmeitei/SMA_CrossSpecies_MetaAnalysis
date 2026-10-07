# Versioned external resources

The analysis does not infer human–mouse equivalence from gene symbols. Before the
cross-species stage, place a versioned Ensembl Compara export at:

`resources/ensembl_human_mouse_one_to_one.tsv`

Required columns are:

`human_ensembl_gene_id`, `human_symbol`, `mouse_ensembl_gene_id`,
`mouse_symbol`, `orthology_type`, `ensembl_release`.

Only rows with `orthology_type == ortholog_one2one` are admitted to the primary
analysis. The helper `R/fetch_orthologs.R` creates this file and records the Ensembl
release.

Network analysis is optional. It requires the complete STRING protein links and
protein-to-gene mapping files. API chunking is intentionally prohibited because it
drops edges between chunks. Record file versions and checksums in
`data/audit/checksums.sha256`.


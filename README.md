#  Gene regulatory networks define human airway epithelial cell types and their distinct responses to Type I interferon
[![DOI](https://img.shields.io/badge/DOI-blue)](https://doi.org/10.64898/2026.05.09.724010)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Codebases and analysis pipelines supporting the manuscript:

> **Gene regulatory networks define human airway epithelial cell types and their distinct responses to Type I interferon** \
> *Bejjani et al.* (2026) *bioRxiv*. DOI: [https://doi.org/10.64898/2026.05.09.724010](https://doi.org/10.64898/2026.05.09.724010)

---

## Interactive Data & Gene Regulatory Network Visualization

Explore the data generated in this study using our Cytoscape sessions and hosted tracks:

* **[Gene Regulatory Network Visualization](https://github.com/MiraldiLab/airwayGRN/blob/main/GRN_viz.md):** Download sessions for interactive visualization of gene regulatory networks using Cytoscape.
* **[UCSC Genome Browser Track Hub](https://github.com/MiraldiLab/airwayGRN/blob/main/genomeBrowserviz.md):** Visualize steady-state and interferon-responsive accessible chromatin and *in silico* [maxATAC](https://doi.org/10.1371/journal.pcbi.1010863) TF binding site predictions, resolved by cell populations and timepoints.

---

## Codebases & Analysis

The analysis is broken down into modular codebases. Detailed instructions for running the code within each module can be found in their respective directories.

* [**TFBS enrichment in IFN-responsive chromatin**](Codebases/Simulate_TFBS_enrich_in_IIPs.R) - Perform simulation-based TFBS enrichment analysis in IFN-increased chromatin regions and accounting for cell type-specific steady-state accessibility.
* [**Enrichment analyses using Fisher's exact test or GSEA**](Codebases/Simulate_TFBS_enrich_in_IIPs.R) - Perform the various enrichment analyses using Fisher's exact tests of GSEA.

---

## Data Availability

Data generated in this manuscript, including Seurat objects, have been deposited on the Gene Expression Ombinus (GEO) and will be made publicly available following publication using the accession numbers:

* scRNA-seq: [GSE330155](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330155)
* snATAC-seq: [GSE330156](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330156)
* snMultiome-seq: [GSE330157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330157)

---

## Additional Resources
* **[maxATAC GitHub Repository](https://github.com/MiraldiLab/maxATAC):** *in silico* TF binding site predictions using [maxATAC](https://doi.org/10.1371/journal.pcbi.1010863).
* **[Inferelator Github Repository](https://github.com/MiraldiLab/InferelatorJL):** GRN inference was performed using the [Inferelator](https://doi.org/10.1186/gb-2006-7-5-r36).
* **[TF-TF Module Analysis Github Repository](https://github.com/MiraldiLab/infTRN_lassoStARS):** Codebase used for TF-TF module analysis (see [`example_Th17_tfTfModules.m`](https://github.com/MiraldiLab/infTRN_lassoStARS/blob/master/Th17_example/example_Th17_tfTfModules.m).
* **[Out-of-Sample Gene Expression Prediction Github Repository](https://github.com/MiraldiLab/infTRN_lassoStARS):** Codebase used for out-of-sample gene expression prediction to determine model complexity (average number of TF regulators per gene, see [`example_workflow_Th17_r2Pred.m`](https://github.com/MiraldiLab/infTRN_lassoStARS/blob/master/Th17_example/example_workflow_Th17_r2Pred.m).

---

## System Requirements

### Software Dependencies
The analyses in this study were performed using R `[4.2.0]` and `[4.2.2]`.

#  Gene regulatory networks define human airway epithelial cell types and their distinct responses to Type I interferon
[![DOI](https://img.shields.io/badge/PENDING_DOI_LINK-blue)](PENDING)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Codebases and analysis pipelines supporting the manuscript:

> **Gene regulatory networks define human airway epithelial cell types and their distinct responses to Type I interferon** \
> *Bejjani et al.* (2026) *Journal Name*. DOI: [PENDING DOI](pending)

---

## Interactive Data & Gene Regulatory Network Visualization

Explore the data generated in this study using our Cytoscape sessions and hosted tracks:

* **[Gene Regulatory Network Visualization](pending/README.md):** Download sessions for interactive visualization of gene regulatory networks using Cytoscape.
* **[UCSC Genome Browser Track Hub](https://github.com/MiraldiLab/airwayGRN/blob/main/genomeBrowserviz.md):** Visualize steady-state and interferon-responsive accessible chromatin and *in silico* [maxATAC](https://doi.org/10.1371/journal.pcbi.1010863) TF binding site predictions, resolved by cell populations and timepoints.

---

## Codebases & Analysis

The analysis is broken down into modular codebases. Detailed instructions for running the code within each module can be found in their respective directories.

* [**TFBS enrichment in IFN-responsive chromatin**](Path/To/Directory/README.md) - Perform TFBS enrichment analysis in IFN-increased chromatin regions and accounting for steady-state accessibility

---

## System Requirements

### Software Dependencies
The analyses in this study were performed using R `[4.2.0]` and `[4.2.2]`. For a full list of package versions, refer to the resource table provided in the manuscript:

[insert table here]

### Additional Resources
* **[maxATAC GitHub Repository](https://github.com/MiraldiLab/maxATAC):** *in silico* TF binding site predictions using [maxATAC](https://doi.org/10.1371/journal.pcbi.1010863).
* **[Inferelator Github Repository](https://github.com/MiraldiLab/InferelatorJL):** GRN inference was performed using the [Inferelator](https://doi.org/10.1186/gb-2006-7-5-r36).
* **[TF-TF module analysis Github Repository](https://github.com/MiraldiLab/infTRN_lassoStARS):** Codebase used for TF-TF module analysis (see [`example_Th17_tfTfModules.m`](https://github.com/MiraldiLab/infTRN_lassoStARS/blob/master/Th17_example/example_Th17_tfTfModules.m).

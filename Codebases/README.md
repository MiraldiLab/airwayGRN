# Instructions

This directory contains the codebases used to process the analyses in the manuscript.

---

## 📂 Directory Setup

The scripts in this folder use relative paths via the `here` package. Before running any code, create the `/inputs` folder and ensure the directory structure looks like this:

```text
📁 Codebases/
├── 📄 README.md
├── 📁 inputs/        
├── 📁 outputs/        
├── 📜 Celltype_IRG_Identification.R
├── 📜 Celltype_IRP_Identification.R
├── 📜 Enrichment_codebases_Fisher_exact_tests.R
├── 📜 Simulate_TFBS_enrich_in_IIPs.R
└── 📜 helper_differential_gene_expression.R
```
The individuals scripts will create `outputs` if you did not do it yourself. Place any inputs in `inputs` or update the paths in the scripts to the correct files.

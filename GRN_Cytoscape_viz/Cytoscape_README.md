## Cytoscape Regulatory Network Visualization
### Session file overview

Explore the GRNs generated in this study using our Cytoscape sessions:

* **[Steady-state core subGRNs](https://github.com/MiraldiLab/airwayGRN/blob/main/GRN_Cytoscape_viz/Steady-state%20cores.cys):** Steady-state (t=0h) cores for basal, suprabasal, ciliated, deuterosomal, ionocyte, and secretory cells. Basal-suprabasal and deuterosomal-ciliated shared networks are also included. Related to Figures 2 and 3.
* **[TF IRG cluster networks](https://github.com/MiraldiLab/airwayGRN/blob/main/GRN_Cytoscape_viz/Celltype%20IRG%20TF%20regulators.cys):** Full and sub-GRNs for the enrichment of TFs in IRG clusters in basal, ciliated and secretory cells. Related to Figure 6.
* **[Shared IRG GRN](https://github.com/MiraldiLab/airwayGRN/blob/main/GRN_Cytoscape_viz/Shared%20IRG%20subnetwork.cys):** Full and sub-GRN describing the regulation of shared IRGs, both IFN-increased and decreased genes. Related to Figure 6.
* **[Mucin and chemokine GRNs](https://github.com/MiraldiLab/airwayGRN/blob/main/GRN_Cytoscape_viz/Chemokine%20and%20mucin%20networks.cys):** Full and sub-GRNs describing mucin and chemokine regulation. Related to Figure 7.

---

## Instructions
To view the gene regulatory subnetworks, you must first install [Cytoscape](https://cytoscape.org). The subnetworks visualized in the manuscript were generated using `version 3.10.4`.
Once you have Cytoscape installed, download the `.cys` sesssion files in the [Cytoscape networks](GRN_Cytoscape_viz) folder and load to a new session as shown below.

<img width="1000" alt="Cytoscape new session" src="https://github.com/MiraldiLab/airwayGRN/blob/main/Diagrams/load_screen.png" />

### Changing heatmap node cell type color
To change the cell type heatmap color for gene expression (targets) and TFA (TFs), navigate to "Image/Chart 1" in the Node tab and select the first option, as shown below.

<img width="1000" alt="Change heatmap color 1" src="https://github.com/MiraldiLab/airwayGRN/blob/main/Diagrams/Change_heatmap_color1.png" />

Remove the currently selected columns and select the columns of interest. **Note that columns need to be added in reverse order (6h timepoint first)**.

<img width="1000" alt="Change heatmap color 2" src="https://github.com/MiraldiLab/airwayGRN/blob/main/Diagrams/Change_heatmap_color2.png" />

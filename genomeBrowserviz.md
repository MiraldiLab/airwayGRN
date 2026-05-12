## UCSC Genome Browser Visualizations
We developed a track data hub to visualize both maxATAC-derived _in-silico_ and chromatin accessibility data for basal, ciliated and secretory cells, at steady-state and in reponse to IFN beta. see [Raney et al., 2014](https://academic.oup.com/bioinformatics/article/30/7/1003/232409) publication.

### Visualization of Genomic Regions
To load tthe track hub, navigate to [Track Data Hubs](https://genome.ucsc.edu/cgi-bin/hgHubConnect?#unlistedHubs), or to "My Data" > "Track Hubs" on the [UCSC Genome Browser website](https://genome.ucsc.edu/index.html). In the URL form, enter `https://gb.research.cchmc.org/hub/group/HAE/hub.txt` to add the track hub.

<img width="1189" height="544" alt="260201_CD4_T_Cell_UCSC_Genome_Browser_Screenshot" src="https://github.com/user-attachments/assets/efd73016-e3dd-4ac9-850d-123aca16e1a8" />

### Visualization of ATAC signal tracks
To modify viewing parameters for the chromatin accessibility tracks or to show specific subpopulations (e.g., the Th1 resting and active subpopulations), navigate to the HAE_Signal_Track settings page for this track and select the desired subopoulations to visualize, as shown below.

<img width="1482" height="882" alt="260201_CD4_T_Cell_UCSC_Genome_Browser_Settings" src="https://github.com/user-attachments/assets/62c322ae-6bb2-43dd-8147-2f20b08fe2c3" />

Once the changes are submitted, only the desired subpopulations appear in the browser, group auto-scaled according to the shown signal tracks.

<img width="2250" height="441" alt="260201_CD4_T_Cell_UCSC_Genome_Browser_Th1_Screenshot" src="https://github.com/user-attachments/assets/7a7394f1-c465-4529-a4d0-675f63dbbce6" />

### Visualization of ATAC signal tracks
By default, TFBS predictions are collapsed for each cell type and time point. You can visualize all TFs with predicted binding sites in that region for the chosen subpopulation by converting the view from "dense" to "full", as shown below.

<img width="1225" height="605" alt="260201_CD4_T_Cell_UCSC_Genome_Browser_Settings_Two" src="https://github.com/user-attachments/assets/9184843d-9c1c-4db3-8c91-eaeffe032694" />

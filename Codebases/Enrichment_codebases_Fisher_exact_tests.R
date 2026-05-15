# Enrichment_codebases_Fisher_exact_tests.R
#
# Code bases pertaining to the enrichment tests performed, as indicated
# by figure numbers and panel letters. For enrichment of TFBSs in cell type-specific
# IFN-responsive peaks (Fig. 5H-M), refer to Simulate_TFBS_enrich_in_IIPs.R
#
# Bejjani et al. (2026) "Gene regulatory networks define human airway epithelial
# cell types and their distinct responses to type I interferon"
# Author: Anthony Bejjani, Cincinnati Children's Hospital Medical Center

# load necessary packages
library(here)
library(fgsea)

# set output directory
outdir <- here('outputs')

# Common data files used in the enrichment analyses
## Cell type IRGs
list_celltype_irgs <- readRDS(here('inputs','celltype_irgs_p5l1.rds'))

## IRG cluster assignments
irg_clusters <- read.delim(here('inputs','cluster_annot.tsv'),sep='\t', header=T, row.names = 1)

## IRP cluster assignments
irp_clusters <- read.delim(here('inputs','cluster_annot_10kb.tsv'), sep='\t', header=T, row.names = 1)

## GRN
network <- read.delim(here('inputs','combined_sp.tsv'), sep='\t', header=F)

## gene signatures
scrna_signatures <- read.table(here('inputs','df_sig_group_major_log2FC0p58_FDR1_all.txt'),sep='\t',header=T)

## counts matrices
scatac_counts <- read.delim(here('inputs','scatac_IFN_bulk_major_group_combat.txt'), sep='\t', header=T, row.names=1)
scrna_counts <- read.delim(here('inputs','scrna_IFN_bulk_major_group_combat.txt'), sep='\t', header=T, row.names=1)

## GO database
GO_db <- readRDS(here('inputs','GOBP_genesets.rds'))
GO_db <- lapply(GO_db, function(x) intersect(x, rownames(counts)))
GO_db <- GO_db[lapply(GO_db, length) >= 5]

GO_db_irg <- lapply(GO_db, function(x) intersect(x, rownames(irg_clusters)))
GO_db_irg <- GO_db_irg[lapply(GO_db_irg, length) >= 5]

bed_ref_peaks <- read.delim(here('inputs','ref_peaks_TSS_2kb_overlap.bed'), sep='\t', header=F, row.names = NULL)
bed_ref_peaks$peak <- paste(bed_ref_peaks$V5,bed_ref_peaks$V6,bed_ref_peaks$V7, sep='-')
bed_ref_peaks <- bed_ref_peaks[,c('peak','V4')]
colnames(bed_ref_peaks) <- c('peak','gene')

##### Fig. 1H - Enrichment of deuterosomal cell peaks in steady-state cell type signatures #####
# 1. Define the cell types and paths
cell_types <- c("Basal", "Suprabasal", "Ciliated", "Deuterosomal", "Ionocyte", "Secretory")

# 2. Pre-filter the signatures to ensure they only contain signature genes
filtered_sigs <- lapply(signatures, function(sig) intersect(sig, rownames(scrna_counts)))
bg_total <- length(rownames(scrna_counts))

# Initialize a list to hold all the results
res_list <- list()

# 3. Loop over each cell type to load its corresponding peak BED file
for (ct_bed in cell_types) {
  
  # Construct the exact file path for intersection of peaks with TSSs
  bed_file <- here('inputs', 'peaks', paste0(ct_bed,'_0_up_intersect.bed'))
  
  if (file.exists(bed_file)) {
    # Load the bed file. V7 contains the gene name from the -b file
    bed_data <- fread(bed_file, header = FALSE, fill = TRUE)
    proximal_genes <- unique(bed_data$V7)
    
    # Intersect the proximal genes with the expressed background universe
    proximal_genes <- intersect(proximal_genes, all_genes)
  } else {
    warning(paste("Could not find file:", bed_file))
    proximal_genes <- character(0)
  }
  
  # 4. Loop over each cell type signature to run the enrichment test
  for (ct_sig in cell_types) {
    sig_genes <- filtered_sigs[[ct_sig]]
    n_sig <- length(sig_genes) # Denominator for percent overlap
    
    # Calculate overlap
    overlap_genes <- intersect(proximal_genes, sig_genes)
    n11 <- length(overlap_genes)
    
    # Calculate the rest of the contingency table variables
    n10 <- n_sig - n11                             # In signature, but NOT near peaks
    n01 <- length(proximal_genes) - n11            # Near peaks, but NOT in signature
    n00 <- bg_total - n11 - n10 - n01              # In neither
    
    # Run the Fisher Exact Test (testing for over-representation / enrichment)
    cont_mat <- matrix(c(n11, n01, n10, n00), nrow = 2)
    f_test <- fisher.test(cont_mat, alternative = "greater")
    
    # Calculate Percent Overlap (using signature as the denominator)
    pct_overlap <- ifelse(n_sig > 0, (n11 / n_sig) * 100, 0)
    
    # Create a row for the data frame
    tmp_df <- data.frame(
      celltype1_bed = ct_bed,
      celltype2_signature = ct_sig,
      p_value = f_test$p.value,
      odds_ratio = as.numeric(f_test$estimate),
      percent_overlap = pct_overlap,
      stringsAsFactors = FALSE
    )
    
    # Append to our list
    res_list[[length(res_list) + 1]] <- tmp_df
  }
}

# 5. Combine all the results into the final data frame
final_enrichment_df <- do.call(rbind, res_list)

# Add a multiple testing correction column
final_enrichment_df$padj <- p.adjust(final_enrichment_df$p_value, method = "BH")
# save results
write.table(final_enrichment_df, file.path(outdir,'final_enrichment_multiome_peaks_in_signatures.tsv'), sep='\t', col.names = T, row.names = F, quote=F)

##### Fig. 2Ai, 2Fi, 3Ai, 3Fi - Enrichment of TF targets in steady-state gene signatures #####
bckgrnd <- unique(c(scrna_signatures[grepl('_0',scrna_signatures$group_major),'Gene'],network$V2))

df_res_cores_t0 <- NULL
for(tf in unique(network$V1)){
  net_targ_up_tmp <- network[network$V1 %in% tf & network$V5 > 0,'V2']
  net_targ_dn_tmp <- network[network$V1 %in% tf & network$V5 < 0,'V2']
  for(cell_sig in c('Basal','Suprabasal','Ciliated','Deuterosomal','Ionocyte','Secretory')){
    signatures_up <- scrna_signatures[scrna_signatures$group_major==paste0(cell_sig,'_0') & scrna_signatures$MeanLog2FC > 0,'Gene']
    signatures_dn <- scrna_signatures[scrna_signatures$group_major==paste0(cell_sig,'_0') & scrna_signatures$MeanLog2FC < 0,'Gene']
    if(length(net_targ_up_tmp) >= 5){
      a <- length(intersect(net_targ_up_tmp, signatures_up))
      b <- length(net_targ_up_tmp) - a
      c <- length(signatures_up) - a
      d <- length(bckgrnd) - (a + b + c)
      contingency_table <- matrix(c(a, b, c, d), 
                                  nrow = 2, 
                                  byrow = TRUE)
      fisher_results <- fisher.test(contingency_table, alternative = "greater")
      df_res_cores_t0 <- rbind(df_res_cores_t0, c(tf,'up',cell_sig, as.numeric(fisher_results$p.value),as.numeric(fisher_results$estimate)))
    }
    if(length(net_targ_dn_tmp) >= 5){
      a <- length(intersect(net_targ_dn_tmp, signatures_dn))
      b <- length(net_targ_dn_tmp) - a
      c <- length(signatures_dn) - a
      d <- length(bckgrnd) - (a + b + c)
      contingency_table <- matrix(c(a, b, c, d), 
                                  nrow = 2, 
                                  byrow = TRUE)
      fisher_results <- fisher.test(contingency_table, alternative = "greater")
      df_res_cores_t0 <- rbind(df_res_cores_t0, c(tf,'dn',cell_sig, as.numeric(fisher_results$p.value),as.numeric(fisher_results$estimate)))
    }
  }
}
df_res_cores_t0 <- as.data.frame(df_res_cores_t0)
colnames(df_res_cores_t0) <- c('tf','dir','signature','pval','OR')
df_res_cores_t0$padj <- 1
for(tf in unique(df_res_cores_t0$tf)){
  df_res_cores_t0[df_res_cores_t0$tf == tf & df_res_cores_t0$dir =='up','padj'] <- p.adjust(as.numeric(df_res_cores_t0[df_res_cores_t0$tf == tf & df_res_cores_t0$dir =='up','pval']), method = 'BH')
  df_res_cores_t0[df_res_cores_t0$tf == tf & df_res_cores_t0$dir =='dn','padj'] <- p.adjust(as.numeric(df_res_cores_t0[df_res_cores_t0$tf == tf & df_res_cores_t0$dir =='dn','pval']), method = 'BH')
}
write.table(df_res_cores, file.path(outdir,'df_res_cores_t0.tsv'), sep='\t', col.names = T, row.names = F, quote=F)

##### Fig. 2B, 2H, 3B - Enrichment of shared TF module targets in GO biological processes #####
networks_folder <- here('inputs','Top7_Networks_clust7')
network_files <- list.files(networks_folder, pattern = "\\min2Targs_sp.tsv$", full.names = FALSE)

bckgrnd <- unique(network[,2])
bckgrnd <- c(bckgrnd, 'ISGF3','STAT1_STAT1')

df_res <- NULL
for(network_file in network_files){
  print(network_file)
  network <- read.delim(file.path(networks_folder, network_file), header=F, sep='\t')
  network_targs <- unique(network[,2])
  df_res_tmp <- NULL
  for(go_geneset in names(GO_db)){
    overlap <- intersect(GO_db[[go_geneset]], network_targs)
    GO_subset <- intersect(GO_db[[go_geneset]], bckgrnd)
    #perform test
    res <- phyper(length(overlap)-1,
                  length(network_targs),
                  length(bckgrnd)-length(network_targs),
                  length(GO_subset),
                  lower.tail= FALSE)
    df_res_tmp <- rbind(df_res_tmp, c(gsub('min2Targs_sp\\.tsv$','',network_file),go_geneset,as.numeric(res)))
    
  }
  df_res_tmp <- as.data.frame(df_res_tmp)
  df_res_tmp$padj <- p.adjust(df_res_tmp[,3], method = 'BH')
  df_res <- rbind(df_res, df_res_tmp)
}
df_res <- as.data.frame(df_res)
colnames(df_res) <- c('module','geneset','pvalue','padj')
write.table(df_res, file.path(networks_folder,'GSEA_GO_full.tsv'), sep='\t', col.names = T, quote = F, row.names = F)

##### Fig. 3J - Enrichment of IRG subclusters in GO biological processes #####
df_irg_enrich_irg_bg <- NULL

for(irg_clust in unique(irg_clusters$across)){
  print(irg_clust)
  irg_clust_up <- rownames(irg_clusters)[irg_clusters$across == irg_clust & irg_clusters$within %in% c(4,3)]
  irg_clust_dn <- rownames(irg_clusters)[irg_clusters$across == irg_clust & irg_clusters$within %in% c(1,2)]
  for(geneset in names(GO_db_irg)){
    if(length(irg_clust_up) > 5){
      
      res <- phyper(length(intersect(GO_db_irg[[geneset]], irg_clust_up))-1,
                    length(irg_clust_up),
                    length(unique(network$V2))-length(irg_clust_up),
                    length(GO_db_irg[[geneset]]),
                    lower.tail= FALSE)
      df_irg_enrich_irg_bg <- rbind(df_irg_enrich_irg_bg, c(geneset,irg_clust,'up', as.numeric(res),length(intersect(GO_db_irg[[geneset]], irg_clust_up))/length(irg_clust_up)))
      
      if(length(irg_clust_dn) > 5){
        res <- phyper(length(intersect(GO_db_irg[[geneset]], irg_clust_dn))-1,
                      length(irg_clust_dn),
                      length(unique(network$V2))-length(irg_clust_dn),
                      length(GO_db_irg[[geneset]]),
                      lower.tail= FALSE)
        df_irg_enrich_irg_bg <- rbind(df_irg_enrich_irg_bg, c(geneset,irg_clust,'dn', as.numeric(res),length(intersect(GO_db_irg[[geneset]], irg_clust_dn))/length(irg_clust_dn)))
      }
    }
  }
}
df_irg_enrich_irg_bg <- as.data.frame(df_irg_enrich_irg_bg)
colnames(df_irg_enrich_irg_bg) <- c('geneset','irg_clust','dir','pval','overlap')
df_irg_enrich_irg_bg$padj <- 1
for(irg_clust in unique(df_irg_enrich_irg_bg$irg_clust)){
  df_irg_enrich_irg_bg[df_irg_enrich_irg_bg$irg_clust == irg_clust & df_irg_enrich_irg_bg$dir =='up','padj'] <- p.adjust(as.numeric(df_irg_enrich_irg_bg[df_irg_enrich_irg_bg$irg_clust == irg_clust & df_irg_enrich_irg_bg$dir =='up','pval']), method = 'BH')
  df_irg_enrich_irg_bg[df_irg_enrich_irg_bg$irg_clust == irg_clust & df_irg_enrich_irg_bg$dir =='dn','padj'] <- p.adjust(as.numeric(df_irg_enrich_irg_bg[df_irg_enrich_irg_bg$irg_clust == irg_clust & df_irg_enrich_irg_bg$dir =='dn','pval']), method = 'BH')
}

##### Fig. 6 - Enrichment of IRG subclusters in genes +/- 2 kb of cell type signature peaks (asterisks) #####
df_res_irg_sig_peaks <- NULL
for(irg_clust in c(2,3,8,7,4,1,6,5)){
  print(irg_clust)
  for(celltype in c('Basal','Ciliated','Secretory')){
    sig_peaks <- readLines(here('inputs', paste0('peaks_sig_major_group_log2FC1p5_FDR10_',celltype,'_0_up.txt')))
    sig_genes <- unique(bed_ref_peaks[bed_ref_peaks$peak %in% sig_peaks,'gene'])
    
    max_vst <- data.frame(max=rowMaxs(as.matrix(scrna_counts[,grepl(celltype, colnames(scrna_counts))])),
                          row.names = rownames(scrna_counts))
    bckgrnd <- rownames(max_vst)[which(max_vst$max >= min(vst_cutoffs[celltype]))]
    
    sig_genes <- intersect(sig_genes, bckgrnd)
    irgs_tmp <- rownames(irg_clusters)[irg_clusters$across == irg_clust]
    irgs_tmp <- intersect(bckgrnd, irgs_tmp)
    res <- phyper(length(intersect(sig_genes, irgs_tmp))-1,
                  length(irgs_tmp),
                  length(bckgrnd)-length(irgs_tmp),
                  length(sig_genes),
                  lower.tail= FALSE)
    df_res_irg_sig_peaks <- rbind(df_res_irg_sig_peaks, c(irg_clust,celltype, as.numeric(res),length(intersect(sig_genes, irgs_tmp))/length(irgs_tmp)))
  }
}
df_res_irg_sig_peaks <- as.data.frame(df_res_irg_sig_peaks)
colnames(df_res_irg_sig_peaks) <- c('IRG','celltype','pval','overlap')
df_res_irg_sig_peaks$padj <- 1
for(irg_clust in 1:8){
  df_res_irg_sig_peaks[df_res_irg_sig_peaks$IRG==irg_clust,'padj'] <- p.adjust(as.numeric(df_res_irg_sig_peaks[df_res_irg_sig_peaks$IRG==irg_clust,'pval']), method = 'BH')
}

##### Fig. 6 - Enrichment of IRG subclusters in cell type signature genes #####
df_res_irg_sig <- NULL
for(irg_clust in c(2,3,8,7,4,1,6,5)){
  print(irg_clust)
  for(celltype in c('Basal','Basal differentiating','Ciliated','FOXN4','Secretory')){
    sig_genes <- readLines(here('inputs', paste0('peaks_sig_group_major_log2FC0p58_FDR10_',celltype,'_0_up.txt')))
    
    celltype <- ifelse(celltype=='Basal differentiating','Suprabasal',celltype)
    max_vst <- data.frame(max=rowMaxs(as.matrix(scrna_counts[,grepl(celltype, colnames(scrna_counts))])),
                          row.names = rownames(scrna_counts))
    bckgrnd <- rownames(max_vst)[which(max_vst$max >= min(vst_cutoffs[celltype]))]
    
    sig_genes <- intersect(bckgrnd,sig_genes)
    irgs_tmp <- rownames(irg_clusters)[irg_clusters$across == irg_clust &
                                         irg_clusters$within %in% ifelse(irg_cluster_direction[as.character(irg_clust)] == 'up', c(4,3),c(1,2))]
    irgs_tmp <- intersect(bckgrnd, irgs_tmp)
    res <- phyper(length(intersect(sig_genes, irgs_tmp))-1,
                  length(irgs_tmp),
                  length(bckgrnd)-length(irgs_tmp),
                  length(sig_genes),
                  lower.tail= FALSE)
    celltype <- ifelse(celltype=='FOXN4','Deuterosomal',celltype)
    df_res_irg_sig <- rbind(df_res_irg_sig, c(irg_clust,celltype, as.numeric(res),length(intersect(sig_genes, irgs_tmp))/length(irgs_tmp)))
  }
}
df_res_irg_sig <- as.data.frame(df_res_irg_sig)
colnames(df_res_irg_sig) <- c('IRG','celltype','pval','overlap')
df_res_irg_sig$padj <- 1
for(irg_clust in 1:8){
  df_res_irg_sig[df_res_irg_sig$IRG==irg_clust,'padj'] <- p.adjust(as.numeric(df_res_irg_sig[df_res_irg_sig$IRG==irg_clust,'pval']), method = 'BH')
}

##### Fig. 6Bi, 6Di - GSEA of TF targets in cell type-specific IFN response #####
list_celltype_irgs2 <- readRDS(here('inputs', 'all_genes_list_igs.rds'))
irgs <- readRDS(here('inputs', 'celltype_irgs_p5l1.rds'))
names(irgs) <- c('Basal','Suprabasal','Ciliated','Deuterosomal','Ionocyte','Secretory')

# get gene rankings
genesets <- list()
for (celltype in c('Basal','Suprabasal','Ciliated', 'Deuterosomal','Secretory')){
  # get target genes of interest
  irgs[[celltype]] <- cbind(irgs[[celltype]], list_celltype_irgs2[[celltype]][rownames(irgs[[celltype]]), paste0('logFC.time_',celltype,'_6_V_',celltype,'_2')])
  colnames(irgs[[celltype]])[ncol(irgs[[celltype]])] <- paste0('logFC.time_',celltype,'_6_V_',celltype,'_2')
  
  # get early irgs (2v0)
  early_irgs <- list_celltype_irgs2[[celltype]]
  early_irgs <- early_irgs[order(-early_irgs[,2]),]
  early_irgs2 <- early_irgs[,2]
  names(early_irgs2) <- rownames(early_irgs)
  genesets[[paste(celltype, 'early', sep='_')]] <- early_irgs2[!grepl(paste(irgs_remove,collapse='|'), names(early_irgs2))]
  
  # get late irgs (6v2)
  late_irgs <- list_celltype_irgs2[[celltype]]
  late_irgs <- late_irgs[order(-late_irgs[,3]),]
  late_irgs2 <- late_irgs[,3]
  names(late_irgs2) <- rownames(late_irgs)
  genesets[[paste(celltype, 'late', sep='_')]] <- late_irgs2[!grepl(paste(irgs_remove,collapse='|'), names(late_irgs2))]
}

tf_genesets <- list()
for (tf in unique(network[,1])){
  tf_genesets[[paste0(tf,'_up')]] <- unique(network[network$V1==tf& network$V5 > 0,'V2'])
  tf_genesets[[paste0(tf,'_dn')]] <- unique(network[network$V1==tf& network$V5 < 0,'V2'])
}

# Run FGSEA
fgseaRes <- list()
for (celltype in c('Basal','Suprabasal','Ciliated', 'Deuterosomal','Secretory')){
  for (time in c('early','late')){
    fgseaRes_tmp <- fgsea(pathways = tf_genesets, 
                          stats    = genesets[[paste(celltype, time, sep='_')]],
                          minSize  = 5,
                          gseaParam = 1,
                          eps      = 0.0)
    saveRDS(fgseaRes_tmp, file.path(outdir, paste0(celltype, '_',time,'_TF_IRG_enrichment.rds')))
    fgseaRes[[paste(celltype, time, sep='_')]] <- as.data.frame(fgseaRes_tmp)
  }
}
saveRDS(fgseaRes, file.path(outdir, 'fgseaRes_TF_IRG_full.rds'))

##### Fig. S3E - Enrichment of IFN-increased peaks in cell type-specific IRG promoters #####
all_promoter_proximal_peaks <- unique(bed_ref_peaks$peak)

list_celltype_irgs <- list_celltype_irgs[c('Basal','Suprabasal','Ciliated','Secretory')]

df_res_irp_irg_enrich <- NULL
for(irp_clust in names(rownames(irp_clusters))){
  print(irp_clust)
  df_res_tmp <- NULL
  overlap_genes_tmp <- bed_ref_peaks[which(bed_ref_peaks$peak %in% rownames(irp_clusters)[irp_clusters$across == irp_clust]),]
  peaks_not_irps <- setdiff(all_promoter_proximal_peaks, unique(overlap_genes_tmp$peak))
  
  for (irg_clust in names(list_celltype_irgs)){
    # get IRGs for that cluster
    irg_tmp <- list_celltype_irgs[[irg_clust]]
    
    # IRPs in proximal promoter of IRG set
    irp_near_irg <- overlap_genes_tmp[overlap_genes_tmp$gene %in% irg_tmp,] 
    # IRPs not in proximal promoter of IRG set
    irp_not_near_irg <- length(unique(overlap_genes_tmp$peak))-length(unique(irp_near_irg$peak)) 
    # non-IRPs in proximal promoter of IRG set
    not_irp_near_irg <- length(unique(bed_ref_peaks[bed_ref_peaks$peak %in% peaks_not_irps &
                                                bed_ref_peaks$gene %in% irg_tmp,'peak']))
    # non-IRPs not in proximal promoter of IRG set
    not_irp_not_near_irg <- length(unique(bed_ref_peaks[bed_ref_peaks$peak %in% peaks_not_irps &
                                                    bed_ref_peaks$gene %in% setdiff(bed_ref_peaks$gene,irg_tmp),'peak']))
    
    hyper_test <- fisher.test(matrix(c(length(unique(irp_near_irg$peak)), #overlap
                                       length(unique(overlap_genes_tmp$peak))-length(unique(irp_near_irg$peak)), 
                                       length(unique(irp_near_irg$peak))+not_irp_near_irg-length(unique(irp_near_irg$peak)),
                                       length(all_promoter_proximal_peaks)), 2, 2), alternative='greater')
    
    df_res_tmp <- rbind(df_res_tmp, c(irg_clust, irp_clust, hyper_test$p.value, hyper_test$estimate, (length(unique(irp_near_irg$peak))/length(ownames(irp_clusters)[irp_clusters$across == irp_clust]))*100))
  }
  df_res_tmp <- cbind(df_res_tmp, p.adjust(df_res_tmp[,3], method = 'BH'))
  df_res_irp_irg_enrich <- rbind(df_res_irp_irg_enrich,df_res_tmp)
}
df_res_irp_irg_enrich <- as.data.frame(df_res_irp_irg_enrich)
colnames(df_res_irp_irg_enrich) <- c('IRG', 'IRP','p_raw','OR','percent', 'p_adj')
df_res_irp_irg_enrich$log10p <- -log10(as.numeric(df_res_irp_irg_enrich$p_adj))
df_res_irp_irg_enrich$log2OR <- log2(as.numeric(df_res_irp_irg_enrich$OR)+1)
write.table(df_res_irp_irg_enrich, file.path(outdir,'df_res_enrichment_unique_IRGs_in_unique_IRPs_2kb_TSS.tsv'), sep='\t', col.names = T, row.names = F, quote=F)

##### Fig. S4D-F - Enrichment of TFBS proximal to IRG subcluster genes #####
#TF motif/TFBS predictions
fimo <- read.table(here('inputs', 'peaks_scatac_IFN_FIMO_res.tsv'), header = T)
fimo$sequence_name <- gsub(':','-',fimo$sequence_name)
fimo <- fimo[,c('sequence_name','motif_alt_id')]

maxatac <- read.table(here('inputs', 'maxatac_concatenated_intersect.bed'), header = F)
maxatac$sequence_name <- paste(maxatac$V7,maxatac$V8,maxatac$V9, sep='-')
maxatac <- maxatac[maxatac$sequence_name %in% ref_peaks,]
maxatac <- maxatac[,c('sequence_name','V6')]
maxatac <- maxatac[!duplicated(maxatac),]
colnames(maxatac) <- c('sequence_name','motif_alt_id')

#replace fimo motifs with maxatac TFs for TFs with maxatac models
fimo <- fimo[!(fimo$motif_alt_id %in% maxatac$motif_alt_id),]
fimo <- rbind(fimo, maxatac)

# cell type peaks
celltype_ref_peaks <- read.delim(here('inputs', 'ref_peaks_celltype_peaks_overlap.bed'), sep='\t', header=F, row.names = NULL)
celltype_ref_peaks$cellpeak <- paste(celltype_ref_peaks$V1,celltype_ref_peaks$V2,celltype_ref_peaks$V3, sep='-')
celltype_ref_peaks$refpeak <- paste(celltype_ref_peaks$V5,celltype_ref_peaks$V6,celltype_ref_peaks$V7, sep='-')
celltype_ref_peaks <- celltype_ref_peaks[,c('cellpeak','refpeak','V4')]
colnames(celltype_ref_peaks) <- c('cellpeak','refpeak','celltype')

all_tfs <- c(intersect(unique(fimo$motif_alt_id),unique(network$V1)),'ISGF3','STAT1_STAT1')
df_tfbs_irgs <- NULL
for(celltype in c('Basal','Ciliated','Secretory')){
  print(celltype)
  max_vst <- data.frame(max=rowMaxs(as.matrix(scrna_counts[,grepl(celltype, colnames(scrna_counts))])),
                        row.names = rownames(scrna_counts))
  expressed_genes <- c(rownames(max_vst)[which(max_vst$max >= min(vst_cutoffs_t0[celltype]))],'ISGF3','STAT1_STAT1')
  expressed_tfs <- intersect(expressed_genes, all_tfs)
  
  bckgrnd <- unique(celltype_ref_peaks[celltype_ref_peaks$celltype == celltype,'refpeak'])
  bckgrnd <- intersect(bckgrnd, bed_ref_peaks$peak)
for(tf in expressed_tfs){
  # get all t=0 peaks with TF motif/TFBS +/- 2kb of any gene - SET 1
  net_targ_tmp <- network[network$V1 %in% tf ,'V2']
  net_targ_tmp <- unique(bed_ref_peaks[bed_ref_peaks$gene %in% net_targ_tmp, 'peak'])
  net_targ_tmp <- intersect(net_targ_tmp, bckgrnd)
  tfbs_peaks_t0 <- unique(fimo[fimo$motif_alt_id == tf & fimo$sequence_name %in% net_targ_tmp,'sequence_name'])
  
    for(irg_clust in unique(irg_clusters$across)){
      #for(dir in c('up','dn')){
      if(length(tfbs_peaks_t0) > 5){
        # get all peaks +/- 2kb of IRG cluster TSS - SET 2
        irgs_tmp <- rownames(irg_clusters)[irg_clusters$across == irg_clust] #& irg_clusters$within %in% ifelse(dir=='up',c(3,4),c(1,2))]
        irgs_tmp <- unique(bed_ref_peaks[bed_ref_peaks$gene %in% irgs_tmp,'peak'])
        irgs_tmp <- unique(intersect(irgs_tmp, bckgrnd))
        if (length(intersect(tfbs_peaks_t0, irgs_tmp)) != 0){
          contingency_table <- matrix(
            c(length(intersect(tfbs_peaks_t0, irgs_tmp)),      length(setdiff(tfbs_peaks_t0, irgs_tmp)),
              length(setdiff(irgs_tmp,tfbs_peaks_t0)), length(setdiff(bckgrnd, c(tfbs_peaks_t0, irgs_tmp)))),
            nrow = 2,
            byrow = TRUE, 
            dimnames = list(
              "TF_Binding"    = c("Has_TFBS_t0", "No_TFBS_t0"),
              "Peak_Location" = c("In_IRG_Cluster", "Not_In_Cluster")
            )
          )
          fisher_res <- fisher.test(contingency_table, alternative = "greater")
        df_tfbs_irgs <- rbind(df_tfbs_irgs, c(tf,irg_clust,celltype, as.numeric(fisher_res$p.value),length(intersect(tfbs_peaks_t0, irgs_tmp))/length(irgs_tmp)))
        }
      }
    }
  }
}

df_tfbs_irgs <- as.data.frame(df_tfbs_irgs)
colnames(df_tfbs_irgs) <- c('tf','irg_clust','celltype','pval','overlap')
df_tfbs_irgs$padj <- 1
for(celltype in unique(df_tfbs_irgs$celltype)){
  for(irg_clust in unique(df_tfbs_irgs$irg_clust)){
  df_tfbs_irgs[df_tfbs_irgs$celltype==celltype & df_tfbs_irgs$irg_clust == irg_clust ,'padj'] <- p.adjust(as.numeric(df_tfbs_irgs[df_tfbs_irgs$celltype==celltype & df_tfbs_irgs$irg_clust == irg_clust ,'pval']), method = 'BH')
}}
write.table(df_tfbs_irgs, file.path(outdir,'df_tfbs_irgs_cellpeaks_allgenes_net_targs.tsv'), sep='\t', col.names = T, row.names = F, quote=F)

##### Fig. S5 - Enrichment of TF targets in IFN-responsive GO biological processes #####
pathways_oi <- c(
  "GOBP_DEFENSE_RESPONSE",
  "GOBP_INNATE_IMMUNE_RESPONSE",
  "GOBP_RESPONSE_TO_VIRUS",
  "GOBP_REGULATION_OF_RESPONSE_TO_BIOTIC_STIMULUS",
  "GOBP_REGULATION_OF_DEFENSE_RESPONSE",
  "GOBP_NEGATIVE_REGULATION_OF_VIRAL_PROCESS",
  "GOBP_REGULATION_OF_IMMUNE_RESPONSE",
  "GOBP_NEGATIVE_REGULATION_OF_VIRAL_GENOME_REPLICATION",
  "GOBP_REGULATION_OF_VIRAL_PROCESS",
  "GOBP_RESPONSE_TO_CYTOKINE",
  "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
  "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
  
  "GOBP_INFLAMMATORY_RESPONSE",
  "GOBP_DEFENSE_RESPONSE",
  "GOBP_RESPONSE_TO_CYTOKINE",
  "GOBP_EPITHELIAL_CELL_DIFFERENTIATION",
  "GOBP_EPITHELIUM_DEVELOPMENT",
  'GOBP_NEGATIVE_REGULATION_OF_MULTICELLULAR_ORGANISMAL_PROCESS',
  "GOBP_NEGATIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION",
  'GOBP_CELL_CELL_SIGNALING',
  'GOBP_NEGATIVE_REGULATION_OF_LOCOMOTION',
  
  "GOBP_TISSUE_MORPHOGENESIS",
  "GOBP_POSITIVE_REGULATION_OF_CELL_DIFFERENTIATION",
  "GOBP_EPHRIN_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_CELLULAR_COMPONENT_MORPHOGENESIS",
  "GOBP_MORPHOGENESIS_OF_AN_EPITHELIUM",
  
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II",
  "GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX",
  "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
  "GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE",
  
  "GOBP_CELL_CYCLE",
  "GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_PROCESS",
  "GOBP_REGULATION_OF_CENTROSOME_CYCLE",
  "GOBP_CELL_DIVISION",
  "GOBP_CELL_CYCLE_DNA_REPLICATION",
  "GOBP_MICROTUBULE_ORGANIZING_CENTER_ORGANIZATION",
  "GOBP_REGULATION_OF_CENTRIOLE_REPLICATION",
  "GOBP_CENTRIOLE_ASSEMBLY",
  "GOBP_MULTI_CILIATED_EPITHELIAL_CELL_DIFFERENTIATION",
  
  "GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE",
  "GOBP_MONOATOMIC_ION_TRANSMEMBRANE_TRANSPORT"
)

my.cluster <- parallel::makeCluster(
  59, 
  type = "PSOCK",
  outfile = ""
)
registerDoParallel(cl = my.cluster)

df_res_tf_irg_bg <- foreach(
  tf = unique(network$V1), 
  .combine = 'rbind', 
  .inorder=FALSE,
  .export=ls(),
  .packages=c('stringr',
              'tidyverse',
              'matrixStats','reshape2')
) %dopar% {
  df_res_tmp <- NULL
  net_targ_up_tmp <- network[network$V1 %in% tf & network$V5 > 0 & network$V2 %in% rownames(irg_clusters),'V2']
  net_targ_dn_tmp <- network[network$V1 %in% tf & network$V5 < 0 & network$V2 %in% rownames(irg_clusters),'V2']
  for(geneset in pathways_oi){
    if(length(net_targ_up_tmp) >= 5){
      res <- phyper(length(intersect(net_targ_up_tmp, GO_db_irg[[geneset]]))-1,
                    length(GO_db_irg[[geneset]]),
                    length(rownames(counts))-length(GO_db_irg[[geneset]]),
                    length(net_targ_up_tmp),
                    lower.tail= FALSE)
      df_res_tmp <- rbind(df_res_tmp, c(tf,'up',geneset, as.numeric(res),length(intersect(net_targ_up_tmp, GO_db_irg[[geneset]]))/length(GO_db_irg[[geneset]])))
    }
    if(length(net_targ_dn_tmp) >= 5){
      res <- phyper(length(intersect(net_targ_dn_tmp, GO_db_irg[[geneset]]))-1,
                    length(GO_db_irg[[geneset]]),
                    length(rownames(counts))-length(GO_db_irg[[geneset]]),
                    length(net_targ_dn_tmp),
                    lower.tail= FALSE)
      df_res_tmp <- rbind(df_res_tmp, c(tf,'dn',geneset, as.numeric(res),length(intersect(net_targ_dn_tmp, GO_db_irg[[geneset]]))/length(GO_db_irg[[geneset]])))
    }
  }
  return(as.data.frame(df_res_tmp))
}
df_res_tf_irg_bg <- as.data.frame(df_res_tf_irg_bg)
colnames(df_res_tf_irg_bg) <- c('tf','dir','geneset','pval','overlap')
df_res_tf_irg_bg$padj <- 1
for(tf in unique(df_res_tf_irg_bg$tf)){
  df_res_tf_irg_bg[df_res_tf_irg_bg$tf == tf & df_res_tf_irg_bg$dir =='up','padj'] <- p.adjust(as.numeric(df_res_tf_irg_bg[df_res_tf_irg_bg$tf == tf & df_res_tf_irg_bg$dir =='up','pval']), method = 'BH')
  df_res_tf_irg_bg[df_res_tf_irg_bg$tf == tf & df_res_tf_irg_bg$dir =='dn','padj'] <- p.adjust(as.numeric(df_res_tf_irg_bg[df_res_tf_irg_bg$tf == tf & df_res_tf_irg_bg$dir =='dn','pval']), method = 'BH')
}
write.table(df_res_tf_irg_bg, file.path(outdir,'df_res_tf_go_hypergeometric_all_genes_bg_all_oi.tsv'), sep='\t', col.names = T, row.names = F, quote=F)

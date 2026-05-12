# Simulate_TFBS_enrich_in_IIPs.R
# simulation-based TFBS enrichment analysis to identify TFs with over- or under- representation 
# of TFBS in IFN-responsive peaks, while controlling for baseline accessibility in the given cell type
# Bejjani et al. (2026) "Gene regulatory networks define human airway epithelial
# cell types and their distinct responses to type I interferon"
# Author: Anthony Bejjani, Cincinnati Children's Hospital Medical Center

# load necessary packages
library(stringr)
library(tidyverse)
library(GenomicRanges)
library(foreach)
library(doParallel)
library(ggplot2)
library(matrixStats)
library(reshape)
library(reshape2)
library(ggrepel)

######## INPUTS ##########
# set output directory
outdir <- '/data/miraldiNB/anthony/projects/HAE/analysis/240811_IRP_sig_peak_enrichment'
dir.create(outdir, recursive = T, showWarnings = F)

# number of parallel processing cores (must be provided for empirical p-value estimation. Cannot run in serial mode; at least 2 cores need to be used)
n.cores <- 60

# gene expression pseudobulk counts
scrna_counts <- read.delim('/data/miraldiNB/anthony/projects/HAE/analysis/230209_pseudobulk_scrna_major/scrna_IFN_bulk_major_group_combat.txt', sep='\t', header=T, row.names=1)
colnames(scrna_counts) <- gsub('Basal_differentiating','Suprabasal', colnames(scrna_counts))
colnames(scrna_counts) <- gsub('FOXN4','Deuterosomal',colnames(scrna_counts))

# chromatin accessibility pseudobulk counts
scatac_counts <- read.delim('/data/miraldiNB/anthony/projects/HAE/analysis/230327_scatac_pseudobulk_2cutsite/scatac_IFN_bulk_major_group_combat.txt', sep='\t', header=T, row.names=1)
colnames(scatac_counts) <- gsub('Basal_differentiating','Suprabasal', colnames(scatac_counts))

# fimo results for all reference peaks
fimo <- read.table("/data/miraldiNB/anthony/projects/HAE/analysis/230824_TF_enrichment_analysis/peaks_scatac_IFN_minvst_5_FIMO_res.tsv", header = T)
fimo$sequence_name <- gsub(':','-',fimo$sequence_name)
tfs <- unique(fimo$motif_alt_id) # get all available TFs

# celltype peaks (or NULL). Will need to be generated from bed files of peaks if not provided.
file_celltype_peaks <- '/data/miraldiNB/anthony/projects/HAE/analysis/240726_IRP_sig_peak_enrichment/celltype_peaks.rds'

# irp clusters and associated cell types. data frame with rows being IRPs and a column of cluster (1 through 6)
irp_clusters <- read.delim('/data/miraldiNB/anthony/projects/HAE/analysis/240313_fig_manuscript_atac/cluster_annot_10kb.tsv', sep='\t', header=T, row.names = 1)

# for each peak cluster, specify the cell types that are relevant for the analysis. Here, we test the enrichment for each cell type background for each cell type.
cluster_ref <- c('1'='Basal, Suprabasal, Ciliated, Deuterosomal, Secretory',
                 '2'='Basal, Suprabasal, Ciliated, Deuterosomal, Secretory',
                 '3'='Basal, Suprabasal, Ciliated, Deuterosomal, Secretory',
                 '4'='Basal, Suprabasal, Ciliated, Deuterosomal, Secretory',
                 '5'='Basal, Suprabasal, Ciliated, Deuterosomal, Secretory',
                 '6'='Basal, Suprabasal, Ciliated, Deuterosomal, Secretory')

######### RUN ##########
# conver reference peaks into GRanges object
ref_peaks <- GRanges(
  seqnames = sapply(strsplit(rownames(scatac_counts), split='-'),"[[",1),
  ranges = IRanges(start=as.numeric(sapply(strsplit(rownames(scatac_counts), split='-'),"[[",2)), end=as.numeric(sapply(strsplit(rownames(scatac_counts), split='-'),"[[",3))),
)

# 1. get peaks per celltype
if(is.null(celltype_peaks)){
  celltype_peaks <- list()
  for (celltype in c('Basal','Basal_differentiating','Ciliated','Secretory')){
    # for each cell type and timepoint, load the peaks and convert them to GRanges objects
    peaks_0 <- read.delim(paste0('/data/miraldiNB/anthony/projects/HAE/analysis/221114_peaks_celltime/peaks/',celltype,'_0_peaks.bed'), sep='\t', header=F)
    peaks_2 <- read.delim(paste0('/data/miraldiNB/anthony/projects/HAE/analysis/221114_peaks_celltime/peaks/',celltype,'_2_peaks.bed'), sep='\t', header=F)
    peaks_6 <- read.delim(paste0('/data/miraldiNB/anthony/projects/HAE/analysis/221114_peaks_celltime/peaks/',celltype,'_6_peaks.bed'), sep='\t', header=F)
    
    peaks_0 <- GRanges(
      seqnames = peaks_0[,1],
      ranges = IRanges(start=peaks_0[,2], end=peaks_0[,3]),
    )
    peaks_2 <- GRanges(
      seqnames = peaks_2[,1],
      ranges = IRanges(start=peaks_2[,2], end=peaks_2[,3]),
    )
    peaks_6 <- GRanges(
      seqnames = peaks_6[,1],
      ranges = IRanges(start=peaks_6[,2], end=peaks_6[,3]),
    )
    celltype_peaks[[celltype]] <- Reduce(union, list(peaks_0, peaks_2, peaks_6)) # find union of all cell type peaks 
    celltype_peaks[[celltype]] <- findOverlaps(celltype_peaks[[celltype]],ref_peaks) # find reference peaks that overlap cell type peaks
    celltype_peaks[[celltype]] <- rownames(scatac_counts)[unique(as.data.frame(celltype_peaks[[celltype]])[,2])] # add them to list as peak names
  }
  names(celltype_peaks) <- c('Basal','Suprabasal','Ciliated','Secretory')
  saveRDS(celltype_peaks, file.path(outdir,'celltype_peaks.rds'))
} else {
  celltype_peaks <- readRDS(file_celltype_peaks) # load existing list if provided
}

# 2. empirical p-value calculation
bins <- seq(4, 10, by=6/10) # create bins of range 4 to 10 as determined from range of max VST counts distribution
labels <- gsub("(?<!^)(\\d{3})$", ",\\1", bins, perl=T)
rangelabels <- paste(head(labels,-1), tail(labels,-1), sep="-") # create bins

# subset ATAC counts to only include t=0 to be used to find VST distribution of IRPs
scatac_counts_0 <- scatac_counts[,grepl('_0', colnames(scatac_counts))]

# setup cluster for parallel processing
# number of cores must be 1 core less than the available cores
my.cluster <- parallel::makeCluster(
  n.cores-1, 
  type = "PSOCK",
  outfile = ""
)
registerDoParallel(cl = my.cluster)

# run empirical p-value pipeline in parallel
# Runs 200 iterations. Can be run multiple times to increase the number of iterations
# Environment is exported to each itration (.export = ls()) to make the necessary files
# available
df_res_empirical <- foreach(
  nsample = 1:200, 
  .combine = 'rbind', 
  .inorder=FALSE,
  .export=ls(),
  .packages=c('stringr',
              'tidyverse',
              'GenomicRanges',
              'matrixStats')
) %dopar% {
  df_res_tmp <- NULL # dataframe to collect empirical results from each iteration
  
  for (celltype in c('Basal','Suprabasal','Ciliated','Secretory')){
    print(paste(nsample,celltype, sep=': '))
    
    # subset ATAC counts to include only t=0 and cell type associated peaks
    scatac_counts_0_celltype <- scatac_counts_0[celltype_peaks[[celltype]],grepl(celltype, colnames(scatac_counts_0))]
    
    # find the VST bin that each peak belongs to (taking the mean across biological replicates)
    scatac_counts_0_celltype$Bin <- cut(rowMeans(scatac_counts_0_celltype), bins, rangelabels)
    
    # find the IRPs that are relevant to that cell type
    celltype_irps <- rownames(irp_clusters)[which(irp_clusters$new_tmp %in% names(cluster_ref)[grepl(celltype,cluster_ref)])]
    
    # find the number of IRPs in each bin based on VST counts
    irp_bins <- table(scatac_counts_0_celltype[celltype_irps,'Bin'])
    
    # get random subsample from all cell type peaks but matching the distribution of IRPs based on the bins. This will be the "IRP" sample
    celltype_irp_rand_sample <- c()
    for(bin in 1:length(irp_bins)){
      celltype_irp_rand_sample <- c(celltype_irp_rand_sample, sample(rownames(scatac_counts_0_celltype)[which(scatac_counts_0_celltype$Bin == names(irp_bins[bin]))], as.numeric(irp_bins[bin]), replace = F))
    }
    
    # get second random subsample from remaining cell type peaks but matching the distribution of IRPs based on the bins. This will be the background sample
    celltype_bg_rand_sample <- c()
    for(bin in 1:length(irp_bins)){
      celltype_bg_rand_sample <- c(celltype_bg_rand_sample, sample(rownames(scatac_counts_0_celltype)[which(scatac_counts_0_celltype$Bin == names(irp_bins[bin]) & rownames(scatac_counts_0_celltype) %in% setdiff(rownames(scatac_counts_0_celltype), celltype_irp_rand_sample))], as.numeric(irp_bins[bin]), replace = F))
    }
    
    # subset fimo results for each sample
    fimo_irp <- fimo[which(fimo$sequence_name %in% celltype_irp_rand_sample),] 
    fimo_bg <- fimo[which(fimo$sequence_name %in% celltype_bg_rand_sample),]
    
    # get number of bp in the random irp set
    n_bp_in_irp_set <- data.frame(peak = celltype_irp_rand_sample)
    n_bp_in_irp_set <- n_bp_in_irp_set %>% separate(col='peak',sep='-', into=c('chr','start','end'), remove = F, convert = T)
    n_bp_in_irp_set <- sum(n_bp_in_irp_set$end - n_bp_in_irp_set$start)
    
    # get number of bp in the random background set
    n_bp_in_bg_set <- data.frame(peak = celltype_bg_rand_sample)
    n_bp_in_bg_set <- n_bp_in_bg_set %>% separate(col='peak',sep='-', into=c('chr','start','end'), remove = F, convert = T)
    n_bp_in_bg_set <- sum(n_bp_in_bg_set$end - n_bp_in_bg_set$start)
    
    # calculate log2OR for each TF in that cell type
    for (tf in tfs){
      # get motif scan results for TF in the random irp set
      motif_in_irp <- fimo_irp[which(fimo_irp$motif_alt_id == tf),]
      
      # find overlapping motifs for TF in the random irp set
      ir <- GRanges(
        seqnames = motif_in_irp$sequence_name,
        ranges = IRanges(motif_in_irp$start, motif_in_irp$stop),
      )
      ir <- range(ir)
      
      # get number of base pairs of motifs in the random irp set for TF
      n_bp_motif_in_irp_set <- sum(ir@ranges@width)
      
      # get motif scan results for TF in the random background set
      motif_in_bg_set <- fimo_bg[which(fimo_bg$motif_alt_id == tf),]
      
      # find overlapping motifs for TF in the random background set
      ir <- GRanges(
        seqnames = motif_in_bg_set$sequence_name,
        ranges = IRanges(motif_in_bg_set$start, motif_in_bg_set$stop),
      )
      ir <- range(ir)
      
      # get number of base pairs of motifs in the random background set for TF
      n_bp_motif_in_bg_set <- sum(ir@ranges@width)
      
      # find odds ratio between random background and random "IRP" cluster
      OR <- (n_bp_motif_in_irp_set/n_bp_in_irp_set)/(n_bp_motif_in_bg_set/n_bp_in_bg_set)
      
      df_res_tmp <- rbind(df_res_tmp, c(celltype, tf, nsample, OR, (n_bp_motif_in_irp_set/n_bp_in_irp_set)*100, (n_bp_motif_in_bg_set/n_bp_in_bg_set)*100))
    }
  }
  return(as.data.frame(df_res_tmp))
}
df_res_empirical <- as.data.frame(df_res_empirical)
colnames(df_res_empirical) <- c('celltype', 'TF','iteration','OR', 'pct_irp','pct_sig')
df_res_empirical$log2OR <- log2(as.numeric(df_res_empirical$OR))
write.table(df_res_empirical, file.path(outdir,'df_res_empirical_200iter.tsv'), sep='\t', col.names = T, row.names = F, quote=F)

# stop cluster
parallel::stopCluster(cl = my.cluster) 

# 3. find empirical OR cutoff that corresponds to a probability of 0.95, such that
# 95% of the OR distribution is below the cutoff
logOR_cutoffs <- data.frame(matrix(nrow=length(unique(df_res_empirical$TF)), ncol=4,0, dimnames = list(unique(df_res_empirical$TF),unique(df_res_empirical$celltype))))
for(tf in unique(df_res_empirical$TF)){
  for(celltype in unique(df_res_empirical$celltype)){
    logOR_cutoffs[tf, celltype] <- quantile(df_res_empirical[which(df_res_empirical$TF==tf & df_res_empirical$celltype == celltype),'log2OR'], probs=0.95, na.rm=TRUE)
  }}
write.table(logOR_cutoffs, file.path(outdir,'df_res_empirical_cutoffs_200iter.tsv'), sep='\t', col.names = T, row.names = F, quote=F)

# optional: plot mean and standard deviation plots by percent background to assess the subsamples
df_tf_comb <- NULL
for(tf in unique(df_res_empirical$TF)){
  for(celltype in unique(df_res_empirical$celltype)){
    if(nrow(df_res_empirical[which(df_res_empirical$celltype == celltype & df_res_empirical$TF == tf),]) > 0){
      df_tf_comb <- rbind(df_tf_comb, c(tf, celltype, 
                                        mean(df_res_empirical[which(df_res_empirical$celltype == celltype & df_res_empirical$TF == tf),'log2OR'][is.finite(df_res_empirical[which(df_res_empirical$celltype == celltype & df_res_empirical$TF == tf),'log2OR'])]),
                                        sd(df_res_empirical[which(df_res_empirical$celltype == celltype & df_res_empirical$TF == tf),'log2OR'][is.finite(df_res_empirical[which(df_res_empirical$celltype == celltype & df_res_empirical$TF == tf),'log2OR'])]),
                                        mean(as.numeric(df_res_empirical[which(df_res_empirical$celltype == celltype & df_res_empirical$TF == tf),'pct_sig'][is.finite(as.numeric(df_res_empirical[which(df_res_empirical$celltype == celltype & df_res_empirical$TF == tf),'pct_sig']))]))))
    }
  }
}
df_tf_comb <- as.data.frame(df_tf_comb)
colnames(df_tf_comb) <- c('TF', 'Celltype','meanOR','sdOR','meanpct_sig')
df_tf_comb$Celltype <- factor(df_tf_comb$Celltype, levels=c('Basal','Suprabasal','Ciliated','Secretory'))

for (plot_var in c('meanOR','sdOR')){
  pdf(file.path(outdir, paste0(plot_var,'_by_pct_sig_empirical.pdf')), height=8, width=6, compress=F)
  ggplot(df_tf_comb, aes(x=as.numeric(meanpct_sig), y=as.numeric(meanOR)))+
    geom_point(alpha=0.5)+
    geom_hline(yintercept = 0, linetype='dashed')+
    facet_wrap(~Celltype, nrow=4)+
    labs(x='% motifs in background', y=ifelse(plot_var == 'meanOR', 'mean log2(OR)', 'SD log2(OR)')) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x))
  dev.off()
}
pdf(file.path(outdir, 'mean_vs_sd.pdf'), height=8, width=6, compress=F)
ggplot(df_tf_comb, aes(y=as.numeric(meanOR), x=as.numeric(sdOR)))+
  geom_point(alpha=0.5)+
  geom_hline(yintercept = 0, linetype='dashed')+
  facet_wrap(~Celltype, nrow=4)+
  labs(y='mean log2(OR)', x='SD log2(OR)')
dev.off()

# optional: plot individual empirical distribution for a TF of interest by celltype
tf_oi <- 'ISGF3'
pdf(file.path(outdir, paste0(celltype,'_log2OR_distributions_empirical.pdf')), compress=F)
ggplot(df_res_empirical[which(df_res_empirical$TF==tf_oi),], aes(x=log2OR, fill=celltype)) +
  geom_histogram(position="identity", alpha=0.5, bins=50)+
  geom_vline(data=melt(logOR_cutoffs[tf_oi,]), aes(xintercept= value, color=variable),
             linetype="dashed")+
  labs(x="log2(OR)", y = "Count")+
  theme_classic()
dev.off()

# 4. find the observed odds ratios
df_res <- NULL
for(irp_clust in 1:6){
  print(irp_clust)
  
  # get relevant cell types for irp cluster
  celltypes <- gsub(' ','', str_split_1(cluster_ref[irp_clust], ','))
  print(celltypes)
  
  # subset fimo results to peaks in that cluster
  fimo_irp <- fimo[which(fimo$sequence_name %in% rownames(irp_clusters)[which(irp_clusters$new_tmp == irp_clust)]),]
  
  # get number of bp in irp cluster
  n_bp_in_irp_set <- data.frame(peak = rownames(irp_clusters)[which(irp_clusters$new_tmp == irp_clust)])
  n_bp_in_irp_set <- n_bp_in_irp_set %>% separate(col='peak',sep='-', into=c('chr','start','end'), remove = F, convert = T)
  n_bp_in_irp_set <- sum(n_bp_in_irp_set$end - n_bp_in_irp_set$start)
  
  df_res_tmp <- NULL
  for (celltype in celltypes){
    print(celltype)
    scatac_counts_0_celltype <- scatac_counts_0[celltype_peaks[[celltype]],grepl(celltype, colnames(scatac_counts_0))]
    scatac_counts_0_celltype$Bin <- cut(rowMeans(scatac_counts_0_celltype), bins, rangelabels)
    celltype_irps <- rownames(irp_clusters)[which(irp_clusters$new_tmp %in% names(cluster_ref)[grepl(celltype,cluster_ref)])] # find relevent irps
    # remove shared IRPs to get a sample of cell type-associated IRP set
    if(irp_clust %in% c(1,2,5,6)){ celltype_irps <- setdiff(celltype_irps, rownames(irp_clusters)[which(irp_clusters$new_tmp %in% c(3,4))])} 
    irp_bins <- table(scatac_counts_0_celltype[celltype_irps,'Bin'])
    
    # get random background set based on IRP vst distribution
    bg_peakset_celltype <- c()
    for(bin in 1:length(irp_bins)){
      bg_peakset_celltype <- c(bg_peakset_celltype, sample(rownames(scatac_counts_0_celltype)[which(scatac_counts_0_celltype$Bin == names(irp_bins[bin]) & rownames(scatac_counts_0_celltype) %in% setdiff(rownames(scatac_counts_0_celltype), rownames(irp_clusters)))], as.numeric(irp_bins[bin]), replace = F))
    }
    
    # subset fimo results to peaks in background peak set
    fimo_bg_peakset <- fimo[which(fimo$sequence_name %in% unique(bg_peakset_celltype)),]
    
    # get number of bp in background peak set
    n_bp_in_bg_set <- data.frame(peak = unique(bg_peakset_celltype))
    n_bp_in_bg_set <- n_bp_in_bg_set %>% separate(col='peak',sep='-', into=c('chr','start','end'), remove = F, convert = T)
    n_bp_in_bg_set <- sum(n_bp_in_bg_set$end - n_bp_in_bg_set$start)
    
    for (tf in tfs){
      # get motif scan results for TF in IRP set
      motif_in_irp <- fimo_irp[which(fimo_irp$motif_alt_id == tf),]
      
      # find overlapping motifs for TF in IRP set
      ir <- GRanges(
        seqnames = motif_in_irp$sequence_name,
        ranges = IRanges(motif_in_irp$start, motif_in_irp$stop),
      )
      ir <- range(ir)
      
      # get number of base pairs of motifs in IRP set for TF
      n_bp_motif_in_irp_set <- sum(ir@ranges@width)
      
      # get motif scan results for TF in background peak set
      motif_in_bg_set <- fimo_bg_peakset[which(fimo_bg_peakset$motif_alt_id == tf),]
      
      # find overlapping motifs for TF in background peak set
      ir <- GRanges(
        seqnames = motif_in_bg_set$sequence_name,
        ranges = IRanges(motif_in_bg_set$start, motif_in_bg_set$stop),
      )
      ir <- range(ir)
      
      # get number of base pairs of motifs in background peak set for TF
      n_bp_motif_in_bg_set <- sum(ir@ranges@width)
      
      # find observed odds ratio between background peak set and IRP set
      OR <- (n_bp_motif_in_irp_set/n_bp_in_irp_set)/(n_bp_motif_in_bg_set/n_bp_in_bg_set)
      
      df_res_tmp <- rbind(df_res_tmp, c(irp_clust, celltype, tf, OR, n_bp_motif_in_irp_set/n_bp_in_irp_set, n_bp_motif_in_bg_set/n_bp_in_bg_set))
    }
  }
  df_res <- rbind(df_res, df_res_tmp)
}
df_res <- as.data.frame(df_res)
colnames(df_res) <- c('IRP','celltype', 'TF','OR', 'pct_irp','pct_sig')
df_res$log2OR <- log2(as.numeric(df_res$OR))
df_res$sign <- sign(df_res$log2OR)
write.table(df_res, file.path(outdir,'odds_irp_motifs_vs_signatures.tsv'), sep='\t', col.names = T, row.names = F, quote=F)
df_res[which(df_res$log2OR == -Inf),'log2OR'] <- 0

# estimate p-values for empirical distributions (Gaussian estimate)
# p-values are corrected using BH method for each cell type across all tests
df_res$pval <- NA
df_res$padj <- NA
for (celltype in c('Basal','Suprabasal','Ciliated','Secretory')){
  for(tf in unique(df_res[which(df_res$celltype == celltype),'TF'])){
    values_tmp <- df_res_empirical[which(df_res_empirical$celltype == celltype & df_res_empirical$TF == tf),'log2OR']
    pvals_tmp <- (1-pnorm(abs(((df_res[which(df_res$celltype == celltype & df_res$TF == tf),'log2OR']-mean(values_tmp[is.finite(values_tmp)]))/sd(values_tmp[is.finite(values_tmp)]))),
                          mean = 0,
                          sd = 1, lower.tail = TRUE))*2
    df_res[which(df_res$celltype == celltype & df_res$TF == tf),'pval'] <- pvals_tmp
  }
  df_res[which(df_res$celltype == celltype),'padj'] <- p.adjust(df_res[which(df_res$celltype == celltype),'pval'], method = 'BH')
}
df_res$log10padj <- -log10(df_res$padj)

# filter results based on empirical OR cutoffs, keeping all timepoints.
# this is done after p-value estimation and correction
df_res_filtered_empirical <- NULL
for(celltype in c('Basal','Suprabasal','Ciliated','Secretory')){
  for (tf in rownames(logOR_cutoffs)[!is.na(logOR_cutoffs[,celltype])]){
    if(!is.na(logOR_cutoffs[tf,celltype])){
      if(nrow(df_res[which(df_res$celltype == celltype & df_res$TF == tf),]) > 0){
        df_res_filtered_empirical <- rbind(df_res_filtered_empirical, df_res[which(df_res$celltype == celltype & df_res$TF == tf & abs(df_res$log2OR) >= logOR_cutoffs[tf,celltype]),])
      }
    }
  }
}
df_res_filtered_empirical <- as.data.frame(df_res_filtered_empirical)
write.table(df_res_filtered_empirical, file.path(outdir,'df_res_filtered_empirical_200iter_all_timepoints.tsv'), sep='\t', col.names = T, row.names = F, quote=F)

###### END OF DATA GENERATION PORTION ######
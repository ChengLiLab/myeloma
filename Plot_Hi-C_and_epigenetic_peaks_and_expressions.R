# Plot Hi-C and epigenetic peaks and expressions

source("../scripts/plotMultipleData2.R")
# read gene expression data
ALL_exp <- read.table("/lustre/user/liclab/wupz/dosageEffect/cuffdiff_GM12878_MM/20160706_cuffdiff_hg19/genes.fpkm_tracking", header = T)
rownames(ALL_exp) <- ALL_exp[,1]
# find gene position
ALL_exp_gene_locus <- strsplit2(ALL_exp[, "locus"], split = ":|-")
rownames(ALL_exp_gene_locus) <- ALL_exp$gene_id

# 1. read ChIP-seq data of U266, chr16
ChIP_names <- c("H3K27ac", "H3K4me1","H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3")
U266_peaks <- list()
for ( i in 1:6 ) {
  U266_peaks[[i]] <- read.table(file = paste( "/lustre/user/liclab/wupz/dosageEffect/plot_hic_epipeaks_expression/U266_", 
                                                            ChIP_names[i], "_peaks.bed", sep = ""),
                                stringsAsFactors = F)
  colnames(U266_peaks[[i]]) <- c("chrom", "start", "end", "peak_name", "height")
  U266_peaks[[i]][U266_peaks[[i]][, 1] == "chrX", 1] <- "chr23"
  U266_peaks[[i]][U266_peaks[[i]][, 1] == "chrY", 1] <- "chr24"
}
head(U266_peaks[[i]])
tail(U266_peaks[[i]])

# 2. Read ChIP-seq data of GM12878, chr16
GM12878_peaks <- list()
for ( i in 1:6 ) {
  GM12878_peaks[[i]] <- read.table(file = paste( "/lustre/user/liclab/wupz/dosageEffect/plot_hic_epipeaks_expression/GM12878_", 
                                              ChIP_names[i], "_peaks.bed", sep = ""),
                                stringsAsFactors = F)
  colnames(GM12878_peaks[[i]]) <- c("chrom", "start", "end", "peak_name", "height")
  GM12878_peaks[[i]][GM12878_peaks[[i]][, 1] == "chrX", 1] <- "chr23"
  GM12878_peaks[[i]][GM12878_peaks[[i]][, 1] == "chrY", 1] <- "chr24"
}
head(GM12878_peaks[[i]])
tail(GM12878_peaks[[i]])

GM12878_peaks_chr16 <- list( )
for ( i in 1:6 ) {
  GM12878_peaks_chr16[[i]] <- read.table(file = paste( "/lustre/user/liclab/wupz/dosageEffect/plot_hic_epipeaks_expression/GM12878_", 
                                                    ChIP_names[i], "_peaks_chr16.bed", sep = "") )
  colnames(GM12878_peaks_chr16[[i]]) <- c("chrom", "start", "end", "peak_name", "height")
}
range(GM12878_peaks_chr16[[3]][, 5])


# 3. Read ICE corrected Hi-C matrix, U266, GM12878, 200kb
U266_HindIII_ICE_cis_matrix
GM12878_HindIII_ICE_cis_matrix <- list()
for ( i in 1:24 ) {
  GM12878_HindIII_ICE_cis_matrix[[i]] <- list()
  tmp_hic_matrix <- read.table( paste("/lustre/user/liclab/lirf/Project/hic/2014/gm12878_hind3_dCTP_296mb/matrix/cis/normalization/chr", 
                                      i, "_normalized_matrix.txt", sep = "") )
  GM12878_HindIII_ICE_cis_matrix[[i]][[i]] <- as.matrix(tmp_hic_matrix)
}
#i = 24
#test <- read.table( paste("/lustre/user/liclab/lirf/Project/hic/2014/gm12878_hind3_dCTP_296mb/resolution_40k/cis/ice_normalization/chr", 
#                          i, "_40k_normalized_matrix.txt", sep = "") )
#dim(test)
#class(test)
# 4. Read Expression of U266 and GM12878

U266_expressions <- data.frame( chrom = as.vector(ALL_exp_gene_locus[, 1]),
                                start = as.numeric(ALL_exp_gene_locus[, 2]),
                                end = as.numeric(ALL_exp_gene_locus[, 3]),
                                gene_id = as.vector(ALL_exp[, "gene_id"]), 
                                U266_FPKM = as.numeric(ALL_exp[, "U266_FPKM"]),
                                strand = ".", 
                                stringsAsFactors = F
                               )
dim(U266_expressions)
head(U266_expressions)
GM12878_expressions <- data.frame( chrom = as.vector(ALL_exp_gene_locus[, 1]),
                                start = as.numeric(ALL_exp_gene_locus[, 2]),
                                end = as.numeric(ALL_exp_gene_locus[, 3]),
                                gene_id = as.vector(ALL_exp[, "gene_id"]), 
                                U266_FPKM = as.numeric(ALL_exp[, "GM12878_FPKM"]),
                                strand = ".", 
                                stringsAsFactors = F
)

# read cnv data
U266_CNV_40kb <- read.table("/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_40kb/U266-merged.sorted.bam.dedup.bam_ratio.txt",
                            stringsAsFactors = F, header = T)
head(U266_CNV_40kb)
dim(U266_CNV_40kb)
GM12878_CNV_40kb <- read.table("/lustre/user/liclab/wupz/DNA_seq/CNV_final/GM12878_CNV/40kb/ERR091571.bam_ratio.txt",
                            stringsAsFactors = F, header = T)
head(GM12878_CNV_40kb)
dim(GM12878_CNV_40kb)
# plot chr16 of GM12878 and U266
setwd("/lustre/user/liclab/wupz/dosageEffect/plot_hic_epipeaks_expression/")

plotMultipleData(Hic_matrix = raw_HiC_matrix[[4]][[16]][[16]], 
                 Chipseq_list = U266_peaks_chr16, 
                 RNAseq_list = list(U266_expressions[U266_expressions[, "chrom"] == "chr16", ]),
                 ylabs = c(ChIP_names, "Expression"),
                 main = "U266, chr16", filename = "U266_hic_chip_expression_chr16.png"
                )
# 5. plot Genes in A2B and B2A
tmp_A2B_genes_exp <- A2B_genes_exp

for ( i in rownames( A2B_genes_exp[ A2B_genes_exp_filter, ] ) ) {
  tmp_gene_name = i
  tmp_gene_chr = as.numeric(strsplit2(ALL_exp_gene_locus[i, 1], split = "chr")[2])
  tmp_gene_start = as.numeric( ALL_exp_gene_locus[i, 2] )
  tmp_gene_end = as.numeric( ALL_exp_gene_locus[i, 3] )
  tmp_gene_start_id <- floor( (tmp_gene_start - 5000000)/200000 )
  tmp_gene_end_id <- floor( (tmp_gene_end + 5000000)/200000 )
  tmp_u266_hic_matrix <- raw_HiC_matrix[[4]][[tmp_gene_chr]][[tmp_gene_chr]][tmp_gene_start_id:tmp_gene_end_id, tmp_gene_start_id:tmp_gene_end_id]
  tmp_u266_chip <- list()
  for ( j in 1:6) {
    tmp_u266_chip[[j]] <-  U266_peaks[[j]][U266_peaks[[j]][, 1] == ALL_exp_gene_locus[i, 1], ]
    tmp_u266_chip[[j]][, c(2,3)] <- tmp_u266_chip[[j]][, c(2,3)]/200000
  }
  tmp_u266_rna <- U266_expressions[U266_expressions[, "chrom"] == ALL_exp_gene_locus[i, 1], ]
  tmp_u266_rna[, c(2,3)] <- tmp_u266_rna[, c(2,3)]/200000
  
  plotMultipleData(Hic_matrix = tmp_u266_hic_matrix, 
                   Chipseq_list = tmp_u266_chip,
                   RNAseq_list = list( tmp_u266_rna ), 
                   ylabs = c(ChIP_names, "Expression"), 
                   main = paste("U266", tmp_gene_name ), 
                   filename = paste("U266 ", tmp_gene_name, ".png", sep = ""), 
                   genomic_range = c( (tmp_gene_start - 5000000)/200000, (tmp_gene_end + 5000000)/200000 )
                  )
}

for ( i in rownames( A2B_genes_exp[ A2B_genes_exp_filter, ] ) ) {
  tmp_gene_name = i
  tmp_genomic_range = ALL_exp_gene_locus[i, 1]
  tmp_gene_start = as.numeric( ALL_exp_gene_locus[i, 2] )
  tmp_gene_end = as.numeric( ALL_exp_gene_locus[i, 3] )
  tmp_gene_start_id <- floor( (tmp_gene_start - 5000000)/200000 )
  tmp_gene_end_id <- floor( (tmp_gene_end + 5000000)/200000 )
  tmp_gm12878_hic_matrix <- GM12878_HindIII_ICE_cis_matrix[[tmp_gene_chr]][[tmp_gene_chr]][tmp_gene_start_id:tmp_gene_end_id, tmp_gene_start_id:tmp_gene_end_id]
  tmp_gm12878_chip <- list()
  for ( j in 1:6) {
    tmp_gm12878_chip[[j]] <-  GM12878_peaks[[j]][GM12878_peaks[[j]][, 1] == ALL_exp_gene_locus[i, 1], ]
    tmp_gm12878_chip[[j]][, c(2,3)] <- tmp_gm12878_chip[[j]][, c(2,3)]/200000
  }
  tmp_gm12878_rna <- GM12878_expressions[GM12878_expressions[, "chrom"] == ALL_exp_gene_locus[i, 1], ]
  tmp_gm12878_rna[, c(2,3)] <- tmp_gm12878_rna[, c(2,3)]/200000
  
  plotMultipleData(Hic_matrix = tmp_gm12878_hic_matrix, 
                   Chipseq_list = tmp_gm12878_chip,
                   RNAseq_list = list( tmp_gm12878_rna ), 
                   ylabs = c(ChIP_names, "Expression"), 
                   main = paste("GM12878", tmp_gene_name ), 
                   filename = paste("GM12878 ", tmp_gene_name, ".png", sep = ""), 
                   genomic_range = c( (tmp_gene_start - 5000000)/200000, (tmp_gene_end + 5000000)/200000 )
  )
}

for ( i in rownames( B2A_genes_exp[ B2A_genes_exp_filter, ] ) ) {
  tmp_gene_name = i
  tmp_gene_chr = as.numeric(strsplit2(ALL_exp_gene_locus[i, 1], split = "chr")[2])
  tmp_gene_start = as.numeric( ALL_exp_gene_locus[i, 2] )
  tmp_gene_end = as.numeric( ALL_exp_gene_locus[i, 3] )
  tmp_gene_start_id <- floor( (tmp_gene_start - 5000000)/200000 )
  tmp_gene_end_id <- floor( (tmp_gene_end + 5000000)/200000 )
  tmp_u266_hic_matrix <- raw_HiC_matrix[[4]][[tmp_gene_chr]][[tmp_gene_chr]][tmp_gene_start_id:tmp_gene_end_id, tmp_gene_start_id:tmp_gene_end_id]
  tmp_u266_chip <- list()
  for ( j in 1:6) {
    tmp_u266_chip[[j]] <-  U266_peaks[[j]][U266_peaks[[j]][, 1] == ALL_exp_gene_locus[i, 1], ]
    tmp_u266_chip[[j]][, c(2,3)] <- tmp_u266_chip[[j]][, c(2,3)]/200000
  }
  tmp_u266_rna <- U266_expressions[U266_expressions[, "chrom"] == ALL_exp_gene_locus[i, 1], ]
  tmp_u266_rna[, c(2,3)] <- tmp_u266_rna[, c(2,3)]/200000
  
  plotMultipleData(Hic_matrix = tmp_u266_hic_matrix, 
                   Chipseq_list = tmp_u266_chip,
                   RNAseq_list = list( tmp_u266_rna ), 
                   ylabs = c(ChIP_names, "Expression"), 
                   main = paste("U266", tmp_gene_name ), 
                   filename = paste("U266 ", tmp_gene_name, ".png", sep = ""), 
                   genomic_range = c( (tmp_gene_start - 5000000)/200000, (tmp_gene_end + 5000000)/200000 )
  )
}

for ( i in rownames( B2A_genes_exp[ B2A_genes_exp_filter, ] ) ) {
  tmp_gene_name = i
  tmp_gene_chr = as.numeric(strsplit2(ALL_exp_gene_locus[i, 1], split = "chr")[2])
  tmp_gene_start = as.numeric( ALL_exp_gene_locus[i, 2] )
  tmp_gene_end = as.numeric( ALL_exp_gene_locus[i, 3] )
  tmp_gene_start_id <- floor( (tmp_gene_start - 5000000)/200000 )
  tmp_gene_end_id <- floor( (tmp_gene_end + 5000000)/200000 )
  tmp_gm12878_hic_matrix <- GM12878_HindIII_ICE_cis_matrix[[tmp_gene_chr]][[tmp_gene_chr]][tmp_gene_start_id:tmp_gene_end_id, tmp_gene_start_id:tmp_gene_end_id]
  tmp_gm12878_chip <- list()
  for ( j in 1:6) {
    tmp_gm12878_chip[[j]] <-  GM12878_peaks[[j]][GM12878_peaks[[j]][, 1] == ALL_exp_gene_locus[i, 1], ]
    tmp_gm12878_chip[[j]][, c(2,3)] <- tmp_gm12878_chip[[j]][, c(2,3)]/200000
  }
  tmp_gm12878_rna <- GM12878_expressions[GM12878_expressions[, "chrom"] == ALL_exp_gene_locus[i, 1], ]
  tmp_gm12878_rna[, c(2,3)] <- tmp_gm12878_rna[, c(2,3)]/200000
  
  plotMultipleData(Hic_matrix = tmp_gm12878_hic_matrix, 
                   Chipseq_list = tmp_gm12878_chip,
                   RNAseq_list = list( tmp_gm12878_rna ), 
                   ylabs = c(ChIP_names, "Expression"), 
                   main = paste("GM12878", tmp_gene_name ), 
                   filename = paste("GM12878 ", tmp_gene_name, ".png", sep = ""), 
                   genomic_range = c( (tmp_gene_start - 5000000)/200000, (tmp_gene_end + 5000000)/200000 )
  )
}

MM_GWAS_snps <- read.table("MM_GWAS_snps.txt", stringsAsFactors = F)
for ( i in 1:7 ) {
  
}

# test plotMultipleData2
plotMultipleData2(Hic_matrix_list = raw_HiC_matrix[[4]], 
                  Chipseq_list = GM12878_peaks,
                  RNAseq_list = list( GM12878_expressions), 
                  ylabs = c(ChIP_names, "Expression"), 
                  main = paste("GM12878", "test"), 
                  genomic_range = "chr16:69000000-76000000"
                  )

# plot genes
# 1. JAK1
ALL_exp_gene_locus["JAK1",]
ALL_exp["MAP4K3", ]
tmp_chr <- ALL_exp_gene_locus["JAK1",1]
tmp_chr_start <- as.integer( ALL_exp_gene_locus["JAK1",2])
tmp_chr_end <- as.integer( ALL_exp_gene_locus["JAK1",2])


plotMultipleData2(Hic_matrix_list = raw_HiC_matrix[[4]], 
                  Chipseq_list = U266_peaks,
                  RNAseq_list = list( U266_expressions), 
                  ylabs = c(ChIP_names, "Expression"), 
                  main = paste("U266", "JAK1"), 
                  genomic_range = paste(tmp_chr, ":", tmp_chr_start-2000000, "-", tmp_chr_end+2000000, sep = ""),
                  filename = paste("U266_", "JAK1", ".png", sep = "")
)


plotMultipleData2(Hic_matrix_list = GM12878_HindIII_ICE_cis_matrix, 
                  Chipseq_list = GM12878_peaks,
                  RNAseq_list = list( GM12878_expressions), 
                  ylabs = c(ChIP_names, "Expression"), 
                  main = paste("GM12878", "JAK1"), 
                  genomic_range = paste(tmp_chr, ":", tmp_chr_start-2000000, "-", tmp_chr_end+2000000, sep = ""),
                  filename = paste("GM12878_", "JAK1", ".png", sep = "")
)

# plot genes
# 2. CBLB, 200kb
ALL_exp_gene_locus["CBLB",]
ALL_exp["CBLB", ]
tmp_chr <- ALL_exp_gene_locus["CBLB",1]
tmp_chr_start <- as.integer( ALL_exp_gene_locus["CBLB",2])
tmp_chr_end <- as.integer( ALL_exp_gene_locus["CBLB",2])
# plot U266 Hi-C matrix
plotMultipleData2(Hic_matrix_list = raw_HiC_matrix[[4]], 
                  Chipseq_list = U266_peaks,
                  RNAseq_list = list( U266_expressions), 
                  ylabs = c(ChIP_names, "Expression"), 
                  main = paste("U266", "CBLB"), 
                  genomic_range = paste(tmp_chr, ":", tmp_chr_start-2000000, "-", tmp_chr_end+2000000, sep = ""),
                  filename = paste("U266_", "CBLB_test20161205", ".png", sep = "")
)

# plot genes
# 2. CBLB, 40kb
ALL_exp_gene_locus["CBLB",]
ALL_exp["CBLB", ]
tmp_chr <- ALL_exp_gene_locus["CBLB",1]
tmp_chr_start <- as.integer( ALL_exp_gene_locus["CBLB",2])
tmp_chr_end <- as.integer( ALL_exp_gene_locus["CBLB",2])
tmp_hic_matrix <- read.table("/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/ice_normalization/chr3_40k_normalized_matrix.txt", stringsAsFactors = F)
tmp_u266_hic_list <- list()
tmp_u266_hic_list[[tmp_chr]] <- list()
tmp_u266_hic_list[[tmp_chr]][[tmp_chr]] <- tmp_hic_matrix
plotMultipleData2(Hic_matrix_list = tmp_u266_hic_list, 
                  Chipseq_list = U266_peaks,
                  RNAseq_list = list( U266_expressions), 
                  ylabs = c(ChIP_names, "Expression"), 
                  main = paste("U266", "CBLB"), resolution = 40000,
                  genomic_range = paste(tmp_chr, ":", tmp_chr_start-2000000, "-", tmp_chr_end+2000000, sep = ""),
                  filename = paste("U266_", "CBLB_test20161205_40kb", ".png", sep = "")
)

plotMultipleData2(Hic_matrix_list = GM12878_HindIII_ICE_cis_matrix, 
                  Chipseq_list = GM12878_peaks,
                  RNAseq_list = list( GM12878_expressions), 
                  ylabs = c(ChIP_names, "Expression"), 
                  main = paste("GM12878", "CBLB"), 
                  genomic_range = paste(tmp_chr, ":", tmp_chr_start-2000000, "-", tmp_chr_end+2000000, sep = ""),
                  filename = paste("GM12878_", "CBLB", ".png", sep = "")
)


# 3. MAP4K3
ALL_exp_gene_locus["MAP4K3",]
ALL_exp["MAP4K3", ]
tmp_chr <- ALL_exp_gene_locus["MAP4K3",1]
tmp_chr_start <- as.integer( ALL_exp_gene_locus["MAP4K3",2])
tmp_chr_end <- as.integer( ALL_exp_gene_locus["MAP4K3",2])


plotMultipleData2(Hic_matrix_list = raw_HiC_matrix[[4]], 
                  Chipseq_list = U266_peaks,
                  RNAseq_list = list( U266_expressions), 
                  ylabs = c(ChIP_names, "Expression"), 
                  main = paste("U266", "MAP4K3"), 
                  genomic_range = paste(tmp_chr, ":", tmp_chr_start-2000000, "-", tmp_chr_end+2000000, sep = ""),
                  filename = paste("U266_", "MAP4K3", ".png", sep = "")
)


plotMultipleData2(Hic_matrix_list = GM12878_HindIII_ICE_cis_matrix, 
                  Chipseq_list = GM12878_peaks,
                  RNAseq_list = list( GM12878_expressions), 
                  ylabs = c(ChIP_names, "Expression"), 
                  main = paste("GM12878", "MAP4K3"), 
                  genomic_range = paste(tmp_chr, ":", tmp_chr_start-2000000, "-", tmp_chr_end+2000000, sep = ""),
                  filename = paste("GM12878_", "MAP4K3", ".png", sep = "")
)

# 3. plot all genes in enrichment analysis(enrichment_results_of_all_genes_with_compartment_switching_and_diff_exp )
enrichment_results_of_all_genes_with_compartment_switching_and_diff_exp <- read.table(
  file = "../comparments_genes/enrichment_results_of_all_genes_with_compartment_switching_and_diff_exp.txt", 
  sep = "\t", header = T)
gene_name_vectors <- enrichment_results_of_all_genes_with_compartment_switching_and_diff_exp[, 6]
gene_name_vectors <- strsplit2(gene_name_vectors, split = ",\ ")

A2B_KEGG_2016_table <- read.table(file = "A2B_KEGG_2016_table.txt", sep = "\t", header = T)
B2A_KEGG_2016_table <- read.table(file = "B2A_KEGG_2016_table.txt", sep = "\t", header = T)
A2B_KEGG_2016_table_pval <- A2B_KEGG_2016_table[A2B_KEGG_2016_table$P.value<=0.05, ] # [1] 21  7
B2A_KEGG_2016_table_pval <- B2A_KEGG_2016_table[B2A_KEGG_2016_table$P.value<=0.05, ]
dim(B2A_KEGG_2016_table_pval) # [1] 7 7
enriched_ab_switching_genes <- rbind(strsplit2(A2B_KEGG_2016_table_pval$Genes, split = ";"), strsplit2(B2A_KEGG_2016_table_pval$Genes, split = ";") )

# combine all genes
gene_name_vectors_all <- vector()
for (i in 1:dim(enriched_ab_switching_genes)[1] ) {
  gene_name_vectors_all <- union(gene_name_vectors_all, enriched_ab_switching_genes[i, ])
}
gene_name_vectors_all <- setdiff(gene_name_vectors_all, "")


for ( i in gene_name_vectors_all ) {
  print(ALL_exp_gene_locus[i,])
  print(ALL_exp[i, ])
  tmp_chr <- ALL_exp_gene_locus[i,1]
  tmp_chr_start <- as.integer( ALL_exp_gene_locus[i,2])
  tmp_chr_end <- as.integer( ALL_exp_gene_locus[i,3])
  if ( tmp_chr == "chrX" ) {
    tmp_chr <- "chr23"
    tmp_chr_id <- 23
  } else if ( tmp_chr == "chrY" ) {
    tmp_chr <- "chr24"
    tmp_chr_id <- 24
  } else if ( !is.null(tmp_chr) ) {
    tmp_chr_id <- as.integer( strsplit2( x = tmp_chr, split = "chr")[2] )
  }
    
  tmp_u266_hic <- read.table( paste( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/ice_normalization/", 
                                     tmp_chr, "_40k_normalized_matrix.txt", sep = "") )
  tmp_u266_hic <- as.matrix(tmp_u266_hic)
  
  tmp_u266_hic_list <- list()
  tmp_u266_hic_list[[tmp_chr_id]] <- list()
  tmp_u266_hic_list[[tmp_chr_id]][[tmp_chr_id]] <- tmp_u266_hic
  tmp_insulation_score <- read.table(paste( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/TAD_boundary/", 
                                            tmp_chr, "_", tmp_chr, "_40k_normalmatrix.txt.is1000001.ids240001.insulation", sep = ""), header = T)
  tmp_insulation_score_tad <- read.table(paste( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/TAD_boundary/", 
                                            tmp_chr, "_", tmp_chr, "_40k_normalmatrix.txt.is1000001.ids240001.insulation.boundaries", sep = ""), header = T)
  try(plotMultipleData2(Hic_matrix_list = tmp_u266_hic_list, 
                    TAD_insulation_score_list = list(tmp_insulation_score, tmp_insulation_score_tad),
                    freec_CNV = U266_CNV_40kb, ploidy = 2,
                    Chipseq_list = U266_peaks,
                    RNAseq_list = list( U266_expressions), 
                    ylabs = c(ChIP_names, "Expression"), 
                    main = paste("U266", i), resolution = 40000,
                    genomic_range = paste(tmp_chr, ":", tmp_chr_start - 2000000, "-", tmp_chr_end + 2000000, sep = ""),
                    filename = paste(i, "_U266_40kb_flank_2M", "_ChIPseq_RNAseq" ,".png", sep = ""),
                    gene_legend = c(i, tmp_chr_start/2 + tmp_chr_end/2)
                    
  ))
  print(paste("Done with U266 gene", i) )
  # plot GM12878 Hi-C and ChIP-seq peaks
  tmp_gm12878_hic <- read.table( paste( "/lustre/user/liclab/lirf/Project/hic/2014/gm12878_hind3_dCTP_296mb/resolution_40k/cis/ice_normalization/", tmp_chr, "_40k_normalized_matrix.txt", sep = "") )
  tmp_gm12878_hic <- as.matrix(tmp_gm12878_hic)
  tmp_gm12878_hic_list <- list()
  tmp_gm12878_hic_list[[tmp_chr_id]] <- list()
  tmp_gm12878_hic_list[[tmp_chr_id]][[tmp_chr_id]] <- tmp_gm12878_hic
  tmp_insulation_score <- read.table(paste( "/lustre/user/liclab/lirf/Project/hic/2014/gm12878_hind3_dCTP_296mb/resolution_40k/cis/TAD_boundary/", 
                                            tmp_chr, "_", tmp_chr, "_40k_normalmatrix.txt.is1000001.ids240001.insulation", sep = ""), header = T)
  tmp_insulation_score_tad <- read.table(paste( "/lustre/user/liclab/lirf/Project/hic/2014/gm12878_hind3_dCTP_296mb/resolution_40k/cis/TAD_boundary/", 
                                                tmp_chr, "_", tmp_chr, "_40k_normalmatrix.txt.is1000001.ids240001.insulation.boundaries", sep = ""), header = T)
  try(plotMultipleData2(Hic_matrix_list = tmp_gm12878_hic_list, 
                    TAD_insulation_score_list = list(tmp_insulation_score, tmp_insulation_score_tad),
                    freec_CNV = GM12878_CNV_40kb, ploidy = 2,
                    Chipseq_list = GM12878_peaks,
                    RNAseq_list = list( GM12878_expressions), 
                    ylabs = c(ChIP_names, "Expression"), 
                    main = paste("GM12878", i), resolution = 40000,
                    genomic_range = paste(tmp_chr, ":", tmp_chr_start - 2000000, "-", tmp_chr_end + 2000000, sep = ""),
                    filename = paste(i, "_GM12878_40kb_flank_2M","_ChIPseq_RNAseq" ,".png", sep = ""),
                    gene_legend = c(i, tmp_chr_start/2 + tmp_chr_end/2)
  ))
  print(paste("Done with GM12878 gene", i) )
}


# Plot Hi-C matrix around a list of SNP associated with MM
for ( i in 1:dim(MM_GWAS_snps)[1] ) {
  #print(ALL_exp_gene_locus[i,])
  #print(ALL_exp[i, ])
  print(MM_GWAS_snps[i, 4])
  tmp_chr <- MM_GWAS_snps[i,1]
  tmp_chr_start <- as.integer( MM_GWAS_snps[i,2])
  tmp_chr_end <- as.integer( MM_GWAS_snps[i,3])
  if ( tmp_chr == "chrX" ) {
    tmp_chr <- "chr23"
    tmp_chr_id <- 23
  } else if ( tmp_chr == "chrY" ) {
    tmp_chr <- "chr24"
    tmp_chr_id <- 24
  } else if ( !is.null(tmp_chr) ) {
    tmp_chr_id <- as.integer( strsplit2( x = tmp_chr, split = "chr")[2] )
  }
  
  tmp_u266_hic <- read.table( paste( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/ice_normalization/", tmp_chr, "_40k_normalized_matrix.txt", sep = "") )
  tmp_u266_hic <- as.matrix(tmp_u266_hic)
  
  tmp_u266_hic_list <- list()
  tmp_u266_hic_list[[tmp_chr_id]] <- list()
  tmp_u266_hic_list[[tmp_chr_id]][[tmp_chr_id]] <- tmp_u266_hic
  plotMultipleData2(Hic_matrix_list = tmp_u266_hic_list, 
                    Chipseq_list = U266_peaks,
                    RNAseq_list = list( U266_expressions), 
                    ylabs = c(ChIP_names, "Expression"), 
                    main = paste("U266", MM_GWAS_snps[i, 4]), resolution = 40000,
                    genomic_range = paste(tmp_chr, ":", tmp_chr_start - 6000000, "-", tmp_chr_end + 6000000, sep = ""),
                    filename = paste(MM_GWAS_snps[i, 4], "_U266", "_ChIPseq_RNAseq" ,".png", sep = ""),
                    gene_legend = c(MM_GWAS_snps[i, 4], tmp_chr_start/2 + tmp_chr_end/2)
                    
  )
  # plot GM12878 Hi-C and ChIP-seq peaks
  tmp_gm12878_hic <- read.table( paste( "/lustre/user/liclab/lirf/Project/hic/2014/gm12878_hind3_dCTP_296mb/resolution_40k/cis/ice_normalization/", tmp_chr, "_40k_normalized_matrix.txt", sep = "") )
  tmp_gm12878_hic <- as.matrix(tmp_gm12878_hic)
  tmp_gm12878_hic_list <- list()
  tmp_gm12878_hic_list[[tmp_chr_id]] <- list()
  tmp_gm12878_hic_list[[tmp_chr_id]][[tmp_chr_id]] <- tmp_gm12878_hic
  plotMultipleData2(Hic_matrix_list = tmp_gm12878_hic_list, 
                    Chipseq_list = GM12878_peaks,
                    RNAseq_list = list( GM12878_expressions), 
                    ylabs = c(ChIP_names, "Expression"), 
                    main = paste("GM12878", MM_GWAS_snps[i, 4]), resolution = 40000,
                    genomic_range = paste(tmp_chr, ":", tmp_chr_start - 6000000, "-", tmp_chr_end + 6000000, sep = ""),
                    filename = paste(MM_GWAS_snps[i, 4], "_GM12878", "_ChIPseq_RNAseq" ,".png", sep = ""),
                    gene_legend = c(MM_GWAS_snps[i, 4], tmp_chr_start/2 + tmp_chr_end/2)
  )
  
}


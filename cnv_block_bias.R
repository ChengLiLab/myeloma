# find CNV block and plot the mean hic interaction score (raw and normed)
# 20160112
# wupz

source("/lustre/user/liclab/wupz/dosageEffect/scripts/FindCNVBlock.R")
source("/lustre/user/liclab/wupz/dosageEffect/scripts/FindLocalHicMatrix.R")
source("/lustre/user/liclab/wupz/dosageEffect/scripts/MedianHicScoreInCNVBlock.R")
# 1. Read data
# 1.1 read genome cnv data
RPMI8226_genome_cnv <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_200kb/8226-merged.sorted.bam.dedup.bam_ratio.txt', 
                                  header = T,  stringsAsFactors = F)
# 1.2 read hic data

# 1.3 read normalization data
# 1.3.1 hicnorm data
load( '/lustre/user/liclab/wupz/dosageEffect/preprocessingData/RPMI8226_HindIII_hicnorm_cis_matrix.RData' )
# 1.3.2 hicnorm with cnv data
load( '/lustre/user/liclab/wupz/dosageEffect/preprocessingData/RPMI8226_HindIII_hicnorm_withcnv_cis_matrix.RData' )
# 1.3.3 raw hic matrix
load('/lustre/user/liclab/wupz/dosageEffect/preprocessingData/raw_HiC_matrix.RData' )
# 1.3.4 ICE norm matrix
load('/lustre/user/liclab/wupz/dosageEffect/preprocessingData/RPMI8226_HindIII_norm_ICEnorm.RData')


# 2. Find CNV block
RPMI8226_cnv_block_200kb <- FindCNVBlock(RPMI8226_genome_cnv )
sapply(RPMI8226_cnv_block_200kb, class)
levels( RPMI8226_cnv_block_200kb$chrom)
# Change chromose X Y to 23, 24
RPMI8226_cnv_block_200kb$chrom[ RPMI8226_cnv_block_200kb$chrom == 'X'] <- 23
RPMI8226_cnv_block_200kb$chrom[ RPMI8226_cnv_block_200kb$chrom == 'Y'] <- 24
# plot the distribution of the length of each block
png(filename = 'Distributiono of length of CNV block.png', width = 1024, height = 1024)
par(cex = 3)
hist(RPMI8226_cnv_block_200kb$chromEnd -  RPMI8226_cnv_block_200kb$chromStart, 
     main = 'Histogram of length of CNV block', xlab = "Length of CNV block")
dev.off()
# plot the distribution of the cnv of each block
png(  filename = 'Distributiono of Copy Number of CNV block in RPMI-8226.png', width = 1024, height = 1024)
par(cex = 3)
hist(RPMI8226_cnv_block_200kb$score, 
     main = 'Histogram of Copy Number of CNV block in RPMI-8226', xlab = "Copy Number of CNV Block")
dev.off()

# 3. get mean of local hic matrix for each cnv block
PlotBiasOfCNV <- function( cnv_block, cnv_block_mean,log2scale = T, legend_site = "bottom", file_name = NULL, ...) {
  # delete zeros and NAs
  tmp_non_zeros <- cnv_block_mean !=0 & !is.na(cnv_block_mean )
  # only show the copy number of 1 ~ 6
  tmp_cnv_1_6_ids <- cnv_block[, 5] >=1 &cnv_block[, 5] <= 6
  if( !is.null(file_name) ) {
    png(filename = file_name, width = 1024, height = 1024)
    par( cex = 3, mgp = c(1.5, 0.5, 0), pch = '.', lwd = 3, mar = c(3, 4, 4, 3))
  }
  tmp_y <- cnv_block_mean[tmp_non_zeros & tmp_cnv_1_6_ids] 
  tmp_x <- cnv_block[tmp_non_zeros & tmp_cnv_1_6_ids, 5]
  if( log2scale) {
    boxplot(log2(tmp_y) ~ tmp_x, col = c('blue', 'green', 'red', 'red', 'red', 'red'), 
            xlab = 'Copy Number of CNV Blocks',
            #ylab = 'Mean Hi-C Interaction Score\n log2 Scale',
            ... )
    # add linear regression line
    
    tmp_lm <- summary( lm( log2(tmp_y) ~ tmp_x ) )
  }
  else if ( !log2scale) {
    boxplot( tmp_y ~ tmp_x, col = c('blue', 'green', 'red', 'red', 'red', 'red'), 
            xlab = 'Copy Number of CNV Blocks',
            #ylab = 'Mean Hi-C Interaction Score\n log2 Scale',
            ... )
    # add linear regression line
    tmp_lm <- summary( lm( tmp_y ~ tmp_x ) )
  }
  tmp_slope <- coef(tmp_lm)[2, 1]
  tmp_intercept <- coef(tmp_lm)[1, 1]
  abline( a = tmp_intercept, b = tmp_slope, col = 'black' )
  legend(legend_site, bty = 'n', cex = 0.8,
         legend = paste('Slope of fitted linear regression lines: ', round(coef(tmp_lm)[2, 1], 4), '\nP value: ', round(coef(tmp_lm)[2, 4], 10) ) 
         )
  if( !is.null(file_name) ) {
    dev.off()
  }
}

# tmp_bin_size <- (RPMI8226_cnv_block_200kb[, 3] - RPMI8226_cnv_block_200kb[, 2] )/200000
# 3.1 Raw matrix
RPMI8226_cnv_block_200kb_raw_mean <- MedianHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, raw_HiC_matrix[[1]],  methods = "mean")
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_raw_mean, main = 'Raw Reads Counts\nRPMI-8226, HindIII',)
               file_name = 'RPMI8226 CNV bias Mean Hi-C Interaction Score (Raw Reads Counts) in CNV Blocks.png')

# 3.2 ICEnorm Hi-C matrix
RPMI8226_cnv_block_200kb_icenorm_mean <- MedianHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_norm_ICEnorm, methods = "mean")
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_icenorm_mean, main = 'ICE Normalization Score\nRPMI-8226, HindIII',)
               file_name = 'RPMI8226 CNV bias Mean Hi-C Interaction Score (ICE Normalization Score) in CNV Block.png')

# 3.3 HiCNorm matrix
RPMI8226_cnv_block_200kb_hicnorm_mean <- MedianHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_hicnorm_cis_matrix,  methods = "mean")
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_hicnorm_mean, main = 'HiCNorm Normalization Score\nRPMI-8226, HindIII',)
               file_name = 'RPMI8226 CNV bias Mean Hi-C Interaction Score (HiCNorm Normalization Score) in CNV Block.png')

# 3.4 HiCNorm with CNV data
RPMI8226_cnv_block_200kb_hicnorm_withcnv_mean <- MedianHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_hicnorm_withcnv_cis_matrix)
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_hicnorm_withcnv_mean, main = 'HiCNorm(With CNV data) Normalization Score\nRPMI-8226, HindIII',)
               file_name = 'RPMI8226 CNV bias Mean Hi-C Interaction Score (HiCNorm (With CNV data) Normalization Score) in CNV Block.png')

# 4. Write the cnv block and median hic score data
RPMI8226_cnv_block_200kb <- cbind(RPMI8226_cnv_block_200kb, 
                                  median_raw_hic = RPMI8226_cnv_block_200kb_raw_mean,
                                  median_icenorm = RPMI8226_cnv_block_200kb_icenorm_mean, 
                                  median_hicnorm = RPMI8226_cnv_block_200kb_hicnorm_mean,
                                  median_hicnorm_withcnv = RPMI8226_cnv_block_200kb_hicnorm_withcnv_mean)
write.table(x = RPMI8226_cnv_block_200kb, file = 'RPMI8226_cnv_block_200kb.txt')

# use diag interactions
# 5.1 Raw matrix
RPMI8226_cnv_block_200kb_raw_diag <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, raw_HiC_matrix[[1]] )
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_raw_diag, main = 'Raw Reads Counts\nRPMI-8226, HindIII', 
               ylab = 'Median Hi-C Interaction Score, Diag\nlog2 Scale', 
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (Raw Reads Counts) in CNV Blocks diag.png')

# 5.2 ICEnorm Hi-C matrix
RPMI8226_cnv_block_200kb_icenorm_diag<- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_norm_ICEnorm)
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_icenorm_diag, main = 'ICE Normalization Score\nRPMI-8226, HindIII',
               ylab = 'Median Hi-C Interaction Score, Diag\nlog2 Scale', 
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (ICE Normalization Score) in CNV Block diag.png')

# 5.3 HiCNorm matrix
RPMI8226_cnv_block_200kb_hicnorm_diag <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_hicnorm_cis_matrix)
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_hicnorm_diag, main = 'HiCNorm Normalization Score\nRPMI-8226, HindIII',
               ylab = 'Median Hi-C Interaction Score, Diag\nlog2 Scale', 
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (HiCNorm Normalization Score) in CNV Block diag.png')

# 5.4 HiCNorm with CNV data
RPMI8226_cnv_block_200kb_hicnorm_withcnv_diag <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_hicnorm_withcnv_cis_matrix)
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_hicnorm_withcnv_diag, main = 'HiCNorm + CNV Normalization Score\nRPMI-8226, HindIII',
               ylab = 'Median Hi-C Interaction Score, Diag\nlog2 Scale', 
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (HiCNorm + CNV data Normalization Score) in CNV Block diag.png')


# 6.0 correct the distance of the normed hic matrix to get the expected hic matrix
RPMI8226_HindIII_hicnorm_withcnv_cis_matrix_expected <- list()
RPMI8226_HindIII_hicnorm_cis_matrix_expect <- list()
RPMI8226_HindIII_norm_ICEnorm_expect <- list()
raw_HiC_matrix_RPMI8226_HindIII_expect <- list()
for( i in 1:24) {
  # 1. hicnorm + cnv
  RPMI8226_HindIII_hicnorm_withcnv_cis_matrix_expected[[i]] <- list()
  RPMI8226_HindIII_hicnorm_withcnv_cis_matrix_expected[[i]][[i]] <- expectHic(RPMI8226_HindIII_hicnorm_withcnv_cis_matrix[[i]][[i]], chr = i)
  # 2. hicnorm
  RPMI8226_HindIII_hicnorm_cis_matrix_expect[[i]] <- list()
  RPMI8226_HindIII_hicnorm_cis_matrix_expect[[i]][[i]] <- expectHic(RPMI8226_HindIII_hicnorm_cis_matrix[[i]][[i]], chr = i)
  # 3. icenorm
  RPMI8226_HindIII_norm_ICEnorm_expect[[i]] <- list()
  RPMI8226_HindIII_norm_ICEnorm_expect[[i]][[i]] <- expectHic(RPMI8226_HindIII_norm_ICEnorm[[i]][[i]], chr = i)
  # 4. raw matrix
  raw_HiC_matrix_RPMI8226_HindIII_expect[[i]] <- list()
  raw_HiC_matrix_RPMI8226_HindIII_expect[[i]][[i]] <- expectHic(raw_HiC_matrix[[1]][[i]][[i]], chr = i)
  
}


# 6.1 Raw matrix, diag
RPMI8226_cnv_block_200kb_raw_diag_expect <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, raw_HiC_matrix_RPMI8226_HindIII_expect)
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_raw_diag_expect, log2scale = F, main = 'Expected Raw Reads Counts\nRPMI-8226, HindIII',)
               ylab = 'Median Hi-C Interaction Score, Diag', legend_site = 'top',
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (Expected Raw Reads Counts) in CNV Blocks diag.png')

# 6.2 ICEnorm Hi-C matrix
RPMI8226_cnv_block_200kb_icenorm_diag_expect <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_norm_ICEnorm_expect)
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_icenorm_diag_expect, log2scale = F, main = 'Expected ICE Normalization Score\nRPMI-8226, HindIII',)
               ylab = 'Median Hi-C Interaction Score, Diag', legend_site = 'top',
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (Expected ICE Normalization Score) in CNV Block diag.png')

# 6.3 HiCNorm matrix
RPMI8226_cnv_block_200kb_hicnorm_diag_expect  <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_hicnorm_cis_matrix_expect)
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_hicnorm_diag_expect, main = 'Expected HiCNorm Normalization Score\nRPMI-8226, HindIII', log2scale = F)
               ylab = 'Median Hi-C Interaction Score, Diag', legend_site = 'top',
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (Expected HiCNorm Normalization Score) in CNV Block diag.png')

# 6.4 HiCNorm with CNV data
RPMI8226_cnv_block_200kb_hicnorm_withcnv_diag_expect  <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_hicnorm_withcnv_cis_matrix_expected )
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_hicnorm_withcnv_diag_expect, log2scale = F)
               main = 'Expected HiCNorm + CNV Normalization Score\nRPMI-8226, HindIII',
               ylab = 'Median Hi-C Interaction Score, Diag', legend_site = 'top',
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (Expected HiCNorm + CNV Normalization Score) in CNV Block diag.png')

# 6.5 Raw matrix
RPMI8226_cnv_block_200kb_raw_median <- MedianHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, raw_HiC_matrix[[1]], methods = 'median' )
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_raw_median, log2scale = F, main = 'Expected Raw Reads Counts\nRPMI-8226, HindIII',)
               ylab = 'Median Hi-C Interaction Score', legend_site = 'top', 
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (Expected Raw Reads Counts) in CNV Blocks.png')

# 6.6 ICEnorm Hi-C matrix
RPMI8226_cnv_block_200kb_icenorm_median <- MedianHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_norm_ICEnorm, methods = 'median' )
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_icenorm_median, log2scale = F, main = 'Expected ICE Normalization Score\nRPMI-8226, HindIII', )
               ylab = 'Median Hi-C Interaction Score', legend_site = 'top', ylim =c(0,6), 
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (Expected ICE Normalization Score) in CNV Block.png')

# 6.7 HiCNorm matrix
RPMI8226_cnv_block_200kb_hicnorm_expect  <- MedianHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_hicnorm_cis_matrix, methods = 'median' )
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_hicnorm_expect, main = 'Expected HiCNorm Normalization Score\nRPMI-8226, HindIII', log2scale = F, 
               ylab = 'Median Hi-C Interaction Score', legend_site = 'top',ylim =c(0,6)
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (Expected HiCNorm Normalization Score) in CNV Block.png')

# 6.8 HiCNorm with CNV data
RPMI8226_cnv_block_200kb_hicnorm_withcnv_median  <- MedianHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_hicnorm_withcnv_cis_matrix, methods = 'median' )
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_hicnorm_withcnv_median, log2scale = F)
               main = 'Expected HiCNorm + CNV Normalization Score\nRPMI-8226, HindIII', 
               ylab = 'Median Hi-C Interaction Score', legend_site = 'top',
               file_name = 'RPMI8226 CNV bias Median Hi-C Interaction Score (Expected HiCNorm + CNV Normalization Score) in CNV Block.png')

# 6.9 Raw matrix, diag
RPMI8226_cnv_block_200kb_raw_diag <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, raw_HiC_matrix[[1]])
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_raw_diag, log2scale = T, main = 'Raw Reads Counts\nRPMI-8226, HindIII',
               ylab = 'Median Hi-C Interaction Score\nDiag bins (log2)', legend_site =  "bottomright",
               file_name = 'RPMI8226 CNV bias raw median diag.png')

# 6.10 ICEnorm Hi-C matrix
RPMI8226_cnv_block_200kb_icenorm_diag <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_norm_ICEnorm)
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_icenorm_diag, log2scale = T, main = 'ICE Normalization Score\nRPMI-8226, HindIII',
               ylab = 'Median Hi-C Interaction Score\nDiag bins (log2)', legend_site =  "bottomright",
               file_name = 'RPMI8226 CNV bias ICE median diag.png')

# 6.11 HiCNorm matrix
RPMI8226_cnv_block_200kb_hicnorm_diag <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_hicnorm_cis_matrix)
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_hicnorm_diag, main = 'HiCNorm Normalization Score\nRPMI-8226, HindIII', log2scale = T,
  ylab = 'Median Hi-C Interaction Score\nDiag bins (log2)', legend_site = "bottomright",
  file_name = 'RPMI8226 CNV bias HiCNorm median diag.png')

# 6.12 HiCNorm with CNV data
RPMI8226_cnv_block_200kb_hicnorm_withcnv_diag  <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_hicnorm_withcnv_cis_matrix )
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_hicnorm_withcnv_diag, log2scale = T, 
  main = 'Expected CNV-HiCNorm Normalization Score\nRPMI-8226, HindIII',
  ylab = 'Median Hi-C Interaction Score\nDiag bins (log2)', legend_site = "bottomright",
  file_name = 'RPMI8226 CNV bias CNV-HiCNorm median diag.png')


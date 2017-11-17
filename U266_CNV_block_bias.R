# find CNV block and plot the mean hic interaction score (raw and normed)
# 20160802, U266
# wupz

source("/lustre/user/liclab/wupz/dosageEffect/scripts/FindCNVBlock.R")
source("/lustre/user/liclab/wupz/dosageEffect/scripts/FindLocalHicMatrix.R")
source("/lustre/user/liclab/wupz/dosageEffect/scripts/MedianHicScoreInCNVBlock.R")
# 1. Read data
# 1.1 read genome cnv data
U266_genome_cnv <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_200kb/U266-merged.sorted.bam.dedup.bam_ratio.txt', 
                                  header = T,  stringsAsFactors = F)
# 1.2 read raw hic data
class(raw_HiC_matrix[[4]])
dim(raw_HiC_matrix[[4]][[1]][[1]]) # [1] 1247 1247
# 1.3 read normalization data
# 1.3.1 hicnorm data
U266_HindIII_hicnorm_cis_matrix <- list()
for (i in 1:24) {
  U266_HindIII_hicnorm_cis_matrix[[i]] <- list()
  U266_HindIII_hicnorm_cis_matrix[[i]][[i]] <- read.table(file = paste('/lustre/user/liclab/wupz/dosageEffect/hicNormaliztion/200kb_diagT/U266_HindIII/hicnorm/cis/chr',i, '_normalized_matrix.txt',sep = ""), 
                                                          stringsAsFactors = F)
}
length(U266_HindIII_hicnorm_cis_matrix)
dim(U266_HindIII_hicnorm_cis_matrix[[1]][[1]]) # [1] 1247 1247
class(U266_HindIII_hicnorm_cis_matrix[[1]][[1]])
U266_HindIII_hicnorm_cis_matrix[[1]][[1]][1:10, 1:10]

# 1.3.2 ICE norm matrix
U266_HindIII_ICE_cis_matrix <- list()
for (i in 1:24) {
  U266_HindIII_ICE_cis_matrix[[i]] <- list()
  U266_HindIII_ICE_cis_matrix[[i]][[i]] <- read.table(file = paste('/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/matrix/cis/ice_normalization/chr',i, '_200kb_normalized_matrix.txt',sep = ""), 
                                                          stringsAsFactors = F)
}
length(U266_HindIII_ICE_cis_matrix)
dim(U266_HindIII_ICE_cis_matrix[[1]][[1]]) # [1] 1247 1247
class(U266_HindIII_ICE_cis_matrix[[1]][[1]])
U266_HindIII_ICE_cis_matrix[[1]][[1]][1:10, 1:10]

# 2. Find CNV block
U266_cnv_block_200kb <- FindCNVBlock(U266_genome_cnv )
sapply(U266_cnv_block_200kb, class)
levels( U266_cnv_block_200kb$chrom) # NULL
# Change chromose X Y to 23, 24
U266_cnv_block_200kb$chrom[ U266_cnv_block_200kb$chrom == 'X'] <- 23
U266_cnv_block_200kb$chrom[ U266_cnv_block_200kb$chrom == 'Y'] <- 24
# plot the distribution of the length of each block
png(filename = 'Distributiono of length of CNV block.png', width = 1024, height = 1024)
par(cex = 3)
hist(U266_cnv_block_200kb$chromEnd -  U266_cnv_block_200kb$chromStart, 
     main = 'Histogram of length of CNV block', xlab = "Length of CNV block")
dev.off()
# plot the distribution of the cnv of each block
png(  filename = 'Distributiono of Copy Number of CNV block in RPMI-8226.png', width = 1024, height = 1024)
par(cex = 3)
hist(U266_cnv_block_200kb$score, 
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
         legend = paste('Slope of fitted linear regression line: ', round(coef(tmp_lm)[2, 1], 4), '\nP value: ', round(coef(tmp_lm)[2, 4], 10) ) 
  )
  if( !is.null(file_name) ) {
    dev.off()
  }
}

# tmp_bin_size <- (U266_cnv_block_200kb[, 3] - U266_cnv_block_200kb[, 2] )/200000
# 3.1 Raw matrix
U266_cnv_block_200kb_raw_median_diag <- DiagHicScoreInCNVBlock( U266_cnv_block_200kb, raw_HiC_matrix[[4]])
PlotBiasOfCNV( U266_cnv_block_200kb, U266_cnv_block_200kb_raw_median_diag, main = 'Raw Reads Counts\nU266, HindIII',
               ylab = 'Median Hi-C Interaction Score\nDiagonal bins (log2 scale)', 
               file_name = 'U266 CNV bias Mean Hi-C Interaction Score (Raw Reads Counts) in CNV Blocks.png')

# 3.2 ICEnorm Hi-C matrix
U266_cnv_block_200kb_icenorm_median_diag <- DiagHicScoreInCNVBlock( U266_cnv_block_200kb, U266_HindIII_ICE_cis_matrix)
PlotBiasOfCNV( U266_cnv_block_200kb, U266_cnv_block_200kb_icenorm_median_diag, main = 'ICE Normalization Score\nU266, HindIII',
               ylab = 'Median Hi-C Interaction Score\nDiagonal bins (log2 scale)', 
               file_name = 'U266 CNV bias Mean Hi-C Interaction Score (ICE Normalization Score) in CNV Block.png')

# 3.3 HiCNorm matrix
U266_cnv_block_200kb_hicnorm_median_diag <- DiagHicScoreInCNVBlock( U266_cnv_block_200kb, U266_HindIII_hicnorm_cis_matrix)
PlotBiasOfCNV( U266_cnv_block_200kb, U266_cnv_block_200kb_hicnorm_median_diag, main = 'HiCNorm Normalization Score\nU266, HindIII',
               ylab = 'Median Hi-C Interaction Score\nDiagonal bins (log2 scale)', 
               file_name = 'U266 CNV bias Mean Hi-C Interaction Score (HiCNorm Normalization Score) in CNV Block.png')


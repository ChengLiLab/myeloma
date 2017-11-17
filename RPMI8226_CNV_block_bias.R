# find CNV block and plot the mean hic interaction score (raw and normed)
# 20160802, RPMI8226
# wupz

source("/lustre/user/liclab/wupz/dosageEffect/scripts/FindCNVBlock.R")
source("/lustre/user/liclab/wupz/dosageEffect/scripts/FindLocalHicMatrix.R")
source("/lustre/user/liclab/wupz/dosageEffect/scripts/MedianHicScoreInCNVBlock.R")
# 1. Read data
# 1.1 read genome cnv data
RPMI8226_genome_cnv <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_200kb/8226-merged.sorted.bam.dedup.bam_ratio.txt', 
                              header = T,  stringsAsFactors = F)
# 1.2 read raw hic data
class(raw_HiC_matrix[[4]])
dim(raw_HiC_matrix[[4]][[1]][[1]]) # [1] 1247 1247
# 1.3 read normalization data
# 1.3.1 hicnorm data
load( '/lustre/user/liclab/wupz/dosageEffect/preprocessingData/RPMI8226_HindIII_hicnorm_cis_matrix.RData' )
length(RPMI8226_HindIII_hicnorm_cis_matrix)
dim(RPMI8226_HindIII_hicnorm_cis_matrix[[1]][[1]]) # [1] 1247 1247
class(RPMI8226_HindIII_hicnorm_cis_matrix[[1]][[1]])
RPMI8226_HindIII_hicnorm_cis_matrix[[1]][[1]][1:10, 1:10]

# 1.3.2 ICE norm matrix
# 1.3.4 ICE norm matrix
load('/lustre/user/liclab/wupz/dosageEffect/preprocessingData/RPMI8226_HindIII_norm_ICEnorm.RData')

length(RPMI8226_HindIII_norm_ICEnorm)
dim(RPMI8226_HindIII_norm_ICEnorm[[1]][[1]]) # [1] 1247 1247
class(RPMI8226_HindIII_norm_ICEnorm[[1]][[1]])
RPMI8226_HindIII_norm_ICEnorm[[1]][[1]][1:10, 1:10]

# 2. Find CNV block
RPMI8226_cnv_block_200kb <- FindCNVBlock(RPMI8226_genome_cnv )
sapply(RPMI8226_cnv_block_200kb, class)
levels( RPMI8226_cnv_block_200kb$chrom) # NULL
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
source("../scripts/PlotBiasOfCNV.r")

# tmp_bin_size <- (RPMI8226_cnv_block_200kb[, 3] - RPMI8226_cnv_block_200kb[, 2] )/200000
# 3.1 Raw matrix
RPMI8226_cnv_block_200kb_raw_median_diag <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, raw_HiC_matrix[[1]])
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_raw_median_diag, main = 'Raw Reads Counts\nRPMI-8226, HindIII',
               ylab = 'Median Hi-C Interaction Score\nDiagonal bins (log2 scale)', 
               file_name = 'RPMI8226 CNV bias Mean Hi-C Interaction Score (Raw Reads Counts) in CNV Blocks.png')

# 3.2 ICEnorm Hi-C matrix
RPMI8226_cnv_block_200kb_icenorm_median_diag <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_norm_ICEnorm)
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_icenorm_median_diag, main = 'ICE Normalization Score\nRPMI-8226, HindIII',
               ylab = 'Median Hi-C Interaction Score\nDiagonal bins (log2 scale)', 
               file_name = 'RPMI8226 CNV bias Mean Hi-C Interaction Score (ICE Normalization Score) in CNV Block.png')

# 3.3 HiCNorm matrix
RPMI8226_cnv_block_200kb_hicnorm_median_diag <- DiagHicScoreInCNVBlock( RPMI8226_cnv_block_200kb, RPMI8226_HindIII_hicnorm_cis_matrix)
PlotBiasOfCNV( RPMI8226_cnv_block_200kb, RPMI8226_cnv_block_200kb_hicnorm_median_diag, main = 'HiCNorm Normalization Score\nRPMI-8226, HindIII',
               ylab = 'Median Hi-C Interaction Score\nDiagonal bins (log2 scale)', 
               file_name = 'RPMI8226 CNV bias Mean Hi-C Interaction Score (HiCNorm Normalization Score) in CNV Block.png')
# CNV Block v.s. insulation Score
# 20160113
# wupz

source("../scripts/FindCNVBlock.R")
source("../scripts/LocalIS.R")
# 1.1 read data of insulation score
RPMI8226_is_40kb <- list()
for ( i in 1:23 ) {
  RPMI8226_is_40kb[[i]] <- read.table( paste('/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release10/resolution_40k/cis/TAD_boundary/chr',
                                             i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation', sep = ''), 
                                       header = T, 
                                       stringsAsFactors = F)
}


# 1.2 read data of TAD boundary
RPMI8226_is_40kb_TAD <- list()
for ( i in 1:23 ) {
  RPMI8226_is_40kb_TAD[[i]] <- read.table( paste( '/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release10/resolution_40k/cis/TAD_boundary/chr', 
                                                   i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation.boundaries', sep = '' ), 
                                            header = T, 
                                            stringsAsFactors = F) 
}


# 1.3 read data of cnv
RPMI8226_cnv_40kb <- read.table( '/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_40kb/8226-merged.sorted.bam.dedup.bam_ratio.txt', 
                                       header = T, 
                                       stringsAsFactors = F)
# Change chromose X Y to 23, 24
RPMI8226_cnv_40kb$Chromosome[ RPMI8226_cnv_40kb$Chromosome == 'X'] <- 23
RPMI8226_cnv_40kb$Chromosome[ RPMI8226_cnv_40kb$Chromosome == 'Y'] <- 24

# 2 Analysis of the overlaping between CNV block boundary and TAD boundary
# 2.1 call CNV block at 40kb resolution
RPMI8226_cnv_block_40kb <- FindCNVBlock( freec_ratio = RPMI8226_cnv_40kb)

# 2.2
# for each CNV block site and TAD site, calculate the insulation score around it.
RPMI8226_CNV_block_40kb_InsulationScore <- numeric()
RPMI8226_CNV_block_40kb_TAD_InsulationScore <- numeric()
for ( i in 1:23) {
  tmp_cnv_block_ins <- t( apply(X = RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ], 
                                MARGIN = 1, 
                                FUN = LocalIS, 
                                insulation_score =  RPMI8226_is_40kb[[i]]) )
  RPMI8226_CNV_block_40kb_InsulationScore <- rbind(RPMI8226_CNV_block_40kb_InsulationScore, tmp_cnv_block_ins)
  
  tmp_TAD_ins <- t( apply( X = RPMI8226_is_40kb_TAD[[i]], 
                           MARGIN = 1, 
                           FUN = LocalIS, 
                           insulation_score = RPMI8226_is_40kb[[i]] ) 
  ) 
  RPMI8226_CNV_block_40kb_TAD_InsulationScore <- rbind(RPMI8226_CNV_block_40kb_TAD_InsulationScore, tmp_TAD_ins)
}
# 2.3 random select the site to calculate the insulation score. random number is the same with TAD numbers
random_ins <- numeric()
for (i in 1:23) {
  tmp_random_number <- dim(RPMI8226_is_40kb_TAD[[i]] )[1]
  tmp_random_site <- RPMI8226_is_40kb[[i]][ sample( x = 1:dim(RPMI8226_is_40kb[[i]])[1], size = tmp_random_number), ]
  tmp_random_ins <- t(apply(tmp_random_site, 1, LocalIS, insulation_score = RPMI8226_is_40kb[[i]] ) )
  random_ins <- rbind(random_ins,  tmp_random_ins)
}

RPMI8226_median_ins_cnv_block <- apply(RPMI8226_CNV_block_40kb_InsulationScore , 2, function(x) median( as.numeric(x), na.rm = T) )
RPMI8226_median_TAD_ins <- apply(RPMI8226_CNV_block_40kb_TAD_InsulationScore, 2, function(x) median( as.numeric(x), na.rm = T) )
median_random_ins <- apply(random_ins, 2, function(x) median( as.numeric(x), na.rm = T) )

# 2.4 plot the insulation score around cnv block boudary
png( filename = 'RPMI8226_CNV_boundaries_insulation_score.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1, 0.2, 0), lwd = 5, tcl =-0.1, mar = c(3,3,3,2))
plot( RPMI8226_median_TAD_ins , type = 'l', 
      xlab = 'Genomic Regions', 
      ylab = 'Median Insulation Score',
      main = 'Insulation Scores Around \nGenomic Regions, RPMI8226',
      ylim = c(min(RPMI8226_median_TAD_ins, RPMI8226_median_ins_cnv_block, median_random_ins, na.rm = T ), 
               max(RPMI8226_median_TAD_ins, RPMI8226_median_ins_cnv_block, median_random_ins, na.rm = T ) ),
      col = "black",
      xlim = c(0, 60), 
      axes = FALSE)
legend('bottomright', bty = 'n', cex = 0.85,
       legend = c('CNV Breakpoints', 'TAD Boudaries', 'Random Breakpoints'),
       col = c('red', 'black', 'blue'), pch = 19)
axis(side = 1, at = c(1, 14, 26, 38, 51), labels = c('-1000kb', '-480kb','0kb','480kb', '1000kb' ))
axis(side = 2)
lines(  RPMI8226_median_ins_cnv_block , col = "red")
lines(  median_random_ins, col = "blue")
dev.off()


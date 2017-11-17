# CNV Block v.s. insulation Score, U266 cell lines
# 20160113
# wupz
# 1.
# 1.1 read data of insulation score
U266_is_40kb <- list()
for ( i in 1:23 ) {
  U266_is_40kb[[i]] <- read.table( paste('/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/TAD_boundary/chr',
                                             i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation', sep = ''), 
                                       header = T, 
                                       stringsAsFactors = F)
}

# 1.2 read data of TAD boundary
U266_is_40kb_TAD <- list()
for ( i in 1:23 ) {
  U266_is_40kb_TAD[[i]] <- read.table( paste( '/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/TAD_boundary/chr', 
                                                  i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation.boundaries', sep = '' ), 
                                           header = T, 
                                           stringsAsFactors = F) 
}

# 1.3 read data of cnv
U266_cnv_40kb <- read.table( '/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_40kb/U266-merged.sorted.bam.dedup.bam_ratio.txt', 
                                 header = T, 
                                 stringsAsFactors = F)
# Change chromose X Y to 23, 24
U266_cnv_40kb$Chromosome[ U266_cnv_40kb$Chromosome == 'X'] <- 23
U266_cnv_40kb$Chromosome[ U266_cnv_40kb$Chromosome == 'Y'] <- 24

# 2 Analysis of the overlaping between CNV block boundary and TAD boundary
# 2.1 call CNV block at 40kb resolution
U266_cnv_block_40kb <- FindCNVBlock( freec_ratio = U266_cnv_40kb)

# for each CNV block site and TAD site, calculate the insulation score around it.
U266_CNV_block_40kb_InsulationScore <- numeric()
U266_CNV_block_40kb_TAD_InsulationScore <- numeric()
for ( i in 1:23) {
  tmp_cnv_block_ins <- Sites_IS( Sites = U266_cnv_block_40kb[U266_cnv_block_40kb[,1] == i, ], insulation_score_list = U266_is_40kb)
  U266_CNV_block_40kb_InsulationScore <- rbind(U266_CNV_block_40kb_InsulationScore, tmp_cnv_block_ins)
  tmp_TAD_ins <- Sites_IS( Sites = U266_is_40kb_TAD[[i]], insulation_score_list = U266_is_40kb) 
  U266_CNV_block_40kb_TAD_InsulationScore <- rbind(U266_CNV_block_40kb_TAD_InsulationScore, tmp_TAD_ins)
}

# random select the site to calculate the insulation score. random number is the same with TAD numbers
U266_random_ins <- numeric()
for (i in 1:23) {
  tmp_random_number <- dim(U266_is_40kb_TAD[[i]] )[1]
  tmp_random_site <- U266_is_40kb[[i]][ sample( x = 1:dim(U266_is_40kb[[i]])[1], size = tmp_random_number), ]
  tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score_list = U266_is_40kb)
  U266_random_ins <- rbind(U266_random_ins,  tmp_random_ins)
}

U266_median_ins_cnv_block <- apply(U266_CNV_block_40kb_InsulationScore , 2, function(x) median( as.numeric(x), na.rm = T) )
U266_median_TAD_ins <- apply(U266_CNV_block_40kb_TAD_InsulationScore, 2, function(x) median( as.numeric(x), na.rm = T) )
U266_median_random_ins <- apply(U266_random_ins, 2, function(x) median( as.numeric(x), na.rm = T) )

# plot the insulation score around cnv block boudary
png( filename = 'U266_CNV_boundaries_insulation_score.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1, 0.2, 0), lwd = 5, tcl =-0.1, mar = c(3,3,3,2))
plot( U266_median_TAD_ins, type = 'l', 
      xlab = 'Genomic Regions', 
      ylab = 'Median Insulation Score',
      main = 'Insulation Scores Around \nGenomic Regions, U266',
      ylim = c(min(U266_median_ins_cnv_block, U266_median_TAD_ins, U266_median_random_ins, na.rm = T), 
               max(U266_median_ins_cnv_block, U266_median_TAD_ins, U266_median_random_ins, na.rm = T) ),
      col = "black",
      xlim = c(0, 60),
      axes = FALSE)
legend('bottomright', bty = 'n', cex = 0.85,
       legend = c('TAD Boudaries', 'CNV Breakpoints', 'Random Sites'),
       col = c('black', 'red', 'blue'), pch = 19)
axis(side = 1, at = c(1, 14, 26, 38, 51), labels = c('-1000kb', '-480kb','0kb','480kb', '1000kb' ))
axis(side = 2)
lines(  U266_median_ins_cnv_block, col = "red")
lines(  U266_median_random_ins, col = "blue")
dev.off()



# Change the tad insulation score from U266 to RPMI-8226
# for each CNV block site and TAD site, calculate the insulation score around it.
U266_CNV_block_40kb_InsulationScore_RPMI8226 <- numeric()
U266_CNV_block_40kb_TAD_InsulationScore_RPMI8226  <- numeric()
for ( i in 1:23) {
  tmp_cnv_block_ins <- Sites_IS( Sites = U266_cnv_block_40kb[U266_cnv_block_40kb[,1] == i, ], insulation_score_list = RPMI8226_is_40kb)
  U266_CNV_block_40kb_InsulationScore_RPMI8226 <- rbind(U266_CNV_block_40kb_InsulationScore_RPMI8226, tmp_cnv_block_ins)
  tmp_TAD_ins <- Sites_IS( Sites = RPMI8226_is_40kb_TAD[[i]], insulation_score_list = RPMI8226_is_40kb) 
  U266_CNV_block_40kb_TAD_InsulationScore_RPMI8226 <- rbind(U266_CNV_block_40kb_TAD_InsulationScore_RPMI8226, tmp_TAD_ins)
}

# random select the site to calculate the insulation score. random number is the same with TAD numbers
random_ins <- numeric()
for (i in 1:23) {
  tmp_random_number <- dim(RPMI8226_is_40kb_TAD[[i]] )[1]
  tmp_random_site <- RPMI8226_is_40kb[[i]][ sample( x = 1:dim(RPMI8226_is_40kb[[i]])[1], size = tmp_random_number), ]
  tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score_list = RPMI8226_is_40kb)
  random_ins <- rbind(random_ins,  tmp_random_ins)
}

median_ins_cnv_block_RPMI8226 <- apply(U266_CNV_block_40kb_InsulationScore_RPMI8226 , 2, function(x) median( as.numeric(x), na.rm = T) )
median_ins_cnv_block <- apply(U266_CNV_block_40kb_InsulationScore , 2, function(x) median( as.numeric(x), na.rm = T) )
median_TAD_ins <- apply(U266_CNV_block_40kb_TAD_InsulationScore_RPMI8226, 2, function(x) median( as.numeric(x), na.rm = T) )
median_random_ins <- apply(random_ins, 2, function(x) median( as.numeric(x), na.rm = T) )

# plot the insulation score around cnv block boudary
png( filename = 'U266_CNV_boundaries_insulation_socre.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1, 0.2, 0), lwd = 3, tcl =-0.1)
plot( median_ins_cnv_block, type = 'l', 
      xlab = 'CNV Block boundary', 
      ylab = 'Median Insulation Score',
      main = 'Insulation Score around boundaries',
      ylim = c(min(median_ins_cnv_block, median_TAD_ins, median_random_ins, na.rm = T), 
               max(median_ins_cnv_block, median_TAD_ins, median_random_ins, na.rm = T) ),
      col = "red",
      axes = FALSE)
legend('bottomright', bty = 'n', cex = 0.65,
       legend = c('CNV Block Boundaries ', 'TAD Boudaries', 'Random Sites'),
       col = c('red', 'blue', 'grey'), lty = 1)
axis(side = 1, at = c(1, 14, 26, 38, 51), labels = c('-1000kb', '-480kb','0kb','480kb', '1000kb' ))
axis(side = 2)

lines(  median_TAD_ins, col = "blue")
lines(  median_random_ins, col = "grey")

dev.off()


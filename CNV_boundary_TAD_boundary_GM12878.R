# CNV Block (RPMI8226, U266) v.s. insulation Score (GM12878), 
# 20160301
# wupz
# 
# 1.1 read data of insulation score
GM12878_is_10kb <- list()
for ( i in c(1:22, 'X') ) {
  GM12878_is_10kb[[i]] <- read.table( paste('/lustre/user/liclab/lirf/Project/hic/disease/data/GM12878_combined/TAD/chr',
                                         i, '_chr', i, '_10kb_normalmatrix.txt.is1000001.ids600001.insulation', sep = ''), 
                                   header = T, 
                                   stringsAsFactors = F)
}
GM12878_is_10kb[[23]] <- GM12878_is_10kb[['X']]


# 1.2 read data of TAD boundary
GM12878_is_10kb_TAD <- list()
for ( i in  c(1:22, 'X')) {
  GM12878_is_10kb_TAD[[i]] <- read.table( paste( '/lustre/user/liclab/lirf/Project/hic/disease/data/GM12878_combined/TAD/chr', 
                                              i, '_chr', i, '_10kb_normalmatrix.txt.is1000001.ids600001.insulation.boundaries', sep = '' ), 
                                       header = T, 
                                       stringsAsFactors = F) 
}
GM12878_is_10kb_TAD[[23]] <- GM12878_is_10kb_TAD[['X']]

# 1.3 read data of cnv, RPMI8226
RPMI8226_cnv_40kb <- read.table( '/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_200kb/8226-merged.sorted.bam.dedup.bam_ratio.txt', 
                                 header = T, 
                                 stringsAsFactors = F)
# Change chromose X Y to 23, 24
RPMI8226_cnv_40kb$Chromosome[ RPMI8226_cnv_40kb$Chromosome == 'X'] <- 23
RPMI8226_cnv_40kb$Chromosome[ RPMI8226_cnv_40kb$Chromosome == 'Y'] <- 24

# 1.4 read data of cnv, U266
U266_cnv_40kb <- read.table( '/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_200kb/U266-merged.sorted.bam.dedup.bam_ratio.txt', 
                             header = T, 
                             stringsAsFactors = F)
# Change chromose X Y to 23, 24
U266_cnv_40kb$Chromosome[ U266_cnv_40kb$Chromosome == 'X'] <- 23
U266_cnv_40kb$Chromosome[ U266_cnv_40kb$Chromosome == 'Y'] <- 24

# 2 Correlation analysis of CNV block boundary and TAD boundary
# 2.1 call CNV block at 40kb resolution, RPMI8226
RPMI8226_cnv_block_40kb <- FindCNVBlock( freec_ratio = RPMI8226_cnv_40kb)
# 2.2 call CNV block at 40kb resolution, U266
U266_cnv_block_40kb <- FindCNVBlock( freec_ratio = U266_cnv_40kb)

# 2.3 for each CNV block sites, TAD sites and random sites, calculate the insulation score around it.
RPMI8226_CNV_block_40kb_InsulationScore_GM12878_10kb <- numeric()
U266_CNV_block_40kb_InsulationScore_GM12878_10kb <- numeric()
TAD_InsulationScore_GM12878_10kb <- numeric()
random_insulationScore <- numeric()

for ( i in 1:23 ) {
  # RPMI8226 CNV block boundary
  print(dim( RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ]))
  tmp_cnv_block_ins <- Sites_IS( Sites = RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ], insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
  RPMI8226_CNV_block_40kb_InsulationScore_GM12878_10kb <- rbind(RPMI8226_CNV_block_40kb_InsulationScore_GM12878_10kb, tmp_cnv_block_ins)
  # U266 CNV block boundary
  print(dim( U266_cnv_block_40kb[U266_cnv_block_40kb[,1] == i, ]))
  tmp_cnv_block_ins <- Sites_IS( Sites = U266_cnv_block_40kb[U266_cnv_block_40kb[,1] == i, ], insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
  U266_CNV_block_40kb_InsulationScore_GM12878_10kb <- rbind(U266_CNV_block_40kb_InsulationScore_GM12878_10kb, tmp_cnv_block_ins)
  # GM12878 TAD boundary
  print(dim(GM12878_is_10kb_TAD[[i]]))
  tmp_TAD_ins <- Sites_IS( Sites = GM12878_is_10kb_TAD[[i]], insulation_score_list = GM12878_is_10kb[[i]], flank = 100) 
  TAD_InsulationScore_GM12878_10kb <- rbind(TAD_InsulationScore_GM12878_10kb, tmp_TAD_ins)
  # Random sites
  print(dim(GM12878_is_10kb_TAD[[i]] ))
  tmp_random_number <- dim(GM12878_is_10kb_TAD[[i]] )[1] # for each chromosome random select some sites with the same number of TADs
  tmp_random_site <- GM12878_is_10kb[[i]][ sample( x = 1:dim(GM12878_is_10kb[[i]])[1], size = tmp_random_number), ]
  tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
  random_insulationScore <- rbind(random_insulationScore,  tmp_random_ins)
}

# 2.4 plot the insulation score around cnv block boudary
png( filename = 'GM12878_TAD_boundary_MM_CNV_boundaries_insulation_socre.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1, 0.2, 0), lwd = 3, tcl =-0.1)
median_cnv_block_is_RMPI8226 <- apply(RPMI8226_CNV_block_40kb_InsulationScore_GM12878_10kb , 2, function(x) median( as.numeric(x), na.rm = T) )
median_cnv_block_is_U266 <- apply(U266_CNV_block_40kb_InsulationScore_GM12878_10kb , 2, function(x) median( as.numeric(x), na.rm = T) )
GM12878_median_TAD_ins <- apply(TAD_InsulationScore_GM12878_10kb, 2, function(x) median( as.numeric(x), na.rm = T) )
GM12878_median_random_ins <- apply(random_insulationScore, 2, function(x) median( as.numeric(x), na.rm = T) )
# plot insulation score around the breakpoints
plot( GM12878_median_TAD_ins, type = 'l', 
      xlab = 'Genomic Sites', 
      ylab = 'Median Insulation Score',
      main = 'CNV Block Boundaries and TAD boundaris',
      ylim = c(min(median_cnv_block_is_RMPI8226, median_cnv_block_is_U266, GM12878_median_TAD_ins,GM12878_median_random_ins, na.rm = T), 
               max(median_cnv_block_is_RMPI8226, median_cnv_block_is_U266, GM12878_median_TAD_ins,GM12878_median_random_ins, na.rm = T) ),
      col = "black",
      axes = FALSE)
lines( GM12878_median_random_ins, col = "grey" )
lines(  median_cnv_block_is_RMPI8226, col = "red")
lines(  median_cnv_block_is_U266, col = "purple")
axis(side = 1, at = c(1, 50, 101, 151, 201), labels = c('-1Mb', '-500kb','0kb','500kb', '1Mb' ))
axis(side = 2)
legend('bottomright', bty = 'n', cex = 0.65,
       legend = c('TAD Boudaries', 'Random Sites','RPMI8226 CNV Block Boundaries', 'U266 CNV Block Boundaries'),
       col = c('black', 'grey', 'red', 'purple'), lty = 1)
dev.off()


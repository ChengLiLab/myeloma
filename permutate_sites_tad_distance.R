# A permutation experiment to test whether the cnv block boundaries is significantly overlapped with TAD boundaries 
# First select random sites, plot the frequency of overlapping sites numbers with TAD boundaries
# Then compare with the frequency of CNV sites (RPMI-8226 and U266 ) that overlapped with TAD boundraies
# 20160329
# wupz
# 1.3 read data of cnv
RPMI8226_cnv_40kb <- read.table( '/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_40kb/8226-merged.sorted.bam.dedup.bam_ratio.txt',
                                 header = T,
                                 stringsAsFactors = F)
# Change chromose X Y to 23, 24
RPMI8226_cnv_40kb$Chromosome[ RPMI8226_cnv_40kb$Chromosome == 'X'] <- 23
RPMI8226_cnv_40kb$Chromosome[ RPMI8226_cnv_40kb$Chromosome == 'Y'] <- 24
RPMI8226_cnv_block_40kb <- FindCNVBlock( freec_ratio = RPMI8226_cnv_40kb)
dim(RPMI8226_cnv_block_40kb)


select_random_sites_one_chrom <- function( all_sites, number_of_selected_sites) {
  select_id <- sample( x = 1:dim(all_sites)[1], size = number_of_selected_sites)
  return( all_sites[select_id, ])
}

# read insulation score data
RPMI8226_is_40kb <- list()
for ( i in 1:23 ) {
  RPMI8226_is_40kb[[i]] <- read.table( paste( '/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release10/resolution_40k/cis/TAD_boundary/chr',
                                                  i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation', sep = '' ),
                                           header = T,
                                           stringsAsFactors = F)
}
# read TAD data
RPMI8226_is_40kb_TAD <- list()
for ( i in 1:23 ) {
  RPMI8226_is_40kb_TAD[[i]] <- read.table( paste( '/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release10/resolution_40k/cis/TAD_boundary/chr',
                                                  i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation.boundaries', sep = '' ),
                                           header = T,
                                           stringsAsFactors = F)
}
length(RPMI8226_is_40kb_TAD)
sum(sapply(RPMI8226_is_40kb_TAD, dim)[1, ])


# find cnv block
source("../scripts/FindCNVBlock.R")
RPMI8226_cnv_block_40kb <- FindCNVBlock( freec_ratio = RPMI8226_cnv_40kb)
# 
RPMI8226_cnv_block_40kb <- RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] != 24, ]
RPMI8226_cnv_block_40kb_tad_distance <- numeric()
source("../scripts/nearestTADDistance.R")
library(GenomicRanges)
for ( i in as.numeric(names(table(RPMI8226_cnv_block_40kb[,1]))) ) {
  tmp_cnv_tad_distance <- nearestTADDistance( query_Sites = RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ],  subject_Sites =  RPMI8226_is_40kb_TAD[[i]])
  RPMI8226_cnv_block_40kb_tad_distance <- c(RPMI8226_cnv_block_40kb_tad_distance, tmp_cnv_tad_distance)
}
RPMI8226_cnv_block_40kb_numbers <- dim(RPMI8226_cnv_block_40kb)[1] # [1] 596
sum(RPMI8226_cnv_block_40kb_tad_distance == 0)/RPMI8226_cnv_block_40kb_numbers # [1] 0.1073826       64/596
sum(RPMI8226_cnv_block_40kb_tad_distance <= 80000)/RPMI8226_cnv_block_40kb_numbers # [1] 0.4194631    250/596
RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb_tad_distance == 0, ]

random_tad_distance <- function(insulation_score, TAD_insulation_score, cnv_block_sites) {
  random_tad_distance <- numeric()
  for ( i in 1:23 ) {
    tmp_random_sites <- select_random_sites_one_chrom(all_sites = insulation_score[[i]], number_of_selected_sites = sum(cnv_block_sites[, 1] == i) )
    tmp_random_sites <- cbind(tmp_random_sites[, 1], tmp_random_sites[, 2] -1)
    tmp_random_tad_distance <- nearestTADDistance( query_Sites = tmp_random_sites,  subject_Sites =  TAD_insulation_score[[i]])
    random_tad_distance <- c(random_tad_distance, tmp_random_tad_distance )
  }
  return(random_tad_distance)
}

# test for 1000 times
# test_random_tad_distance_1000: 40 kb
test_random_tad_distance_1000_RPMI8226 <- numeric()
for ( i in 1:1000) {
  if (i %% 50 == 0) { print(i) }
  test_random_tad_distance_1000_RPMI8226 <- cbind( test_random_tad_distance_1000_RPMI8226, 
                                       random_tad_distance( insulation_score = RPMI8226_is_40kb, TAD_insulation_score = RPMI8226_is_40kb_TAD, cnv_block_sites = RPMI8226_cnv_block_40kb) 
  )
}

# figure 3.c, not used since 2010815
png(filename = "../CNV_TAD_boundary/RPMI8226 Random regions overlapped with TAD Boundaries.png", width = 1024, height = 1024)
par( cex = 4, mgp = c(1.5, 0.5, 0), pch = '.', lwd = 4, mar = c(3, 4, 4, 3))
tmp_hist1 <- hist( apply(test_random_tad_distance_1000_RPMI8226, 2, function(x) sum(x<=80000)), breaks = 18, 
      main = "Random regions overlapped\n with TAD Boundaries, RPMI-8226", xlab = "Permutation times",
      ylab = 'Frequency', col = 'deepskyblue1', cex.axis = 0.8, xlim = c(100, 260) )
arrows(x0 = 250, y0 = 150, x1 = 250, y1 = 0, col = 'red', lwd = 5)
tmp_pval <- sum(apply(test_random_tad_distance_1000_RPMI8226, 2, function(x) sum(x<=80000)) >= 250)/1000
text( x = 225, y = 180, cex = 0.6, col = "red",
      labels = paste('CNV block boundaries \noverlapped with \nTAD Boundaries\nP value:', tmp_pval))
dev.off()


# figure 3.c, modified 2010815
png(filename = "../CNV_TAD_boundary/RPMI8226 Distance of CNV boundaries to TAD Boundaries.png", width = 1024, height = 1024)
par( cex = 3, mgp = c(2, 0.5, 0), pch = '.', lwd = 4, mar = c(5, 4, 4, 3))
RPMI8226_cnv_block_40kb_tad_distance_counts <- table(cut(RPMI8226_cnv_block_40kb_tad_distance, breaks = c(0:25*40000, Inf)-1, include.lowest = T, right = F))
test_random_tad_distance_1000_RPMI8226_counts <- apply(X = test_random_tad_distance_1000_RPMI8226, 2, FUN = function(x) table(cut(x, breaks = c(0:25*40000, Inf)-1, include.lowest = T, right = F)))
test_random_tad_distance_1000_RPMI8226_counts_mean <- apply(test_random_tad_distance_1000_RPMI8226_counts, 1, mean)
test_random_tad_distance_1000_RPMI8226_counts_sd <- apply(test_random_tad_distance_1000_RPMI8226_counts, 1, sd)
barplot_data <- rbind(RPMI8226_cnv_block_40kb_tad_distance_counts, 
                      test_random_tad_distance_1000_RPMI8226_counts_mean)
fisher_test_table <- matrix(c(sum(RPMI8226_cnv_block_40kb_tad_distance_counts[1:3] ), 
                              sum(RPMI8226_cnv_block_40kb_tad_distance_counts[4:26] ), 
                              sum(test_random_tad_distance_1000_RPMI8226_counts_mean[1:3] ),
                              sum(test_random_tad_distance_1000_RPMI8226_counts_mean[4:26])), 
                            nrow = 2, byrow = T)
fisher.test(fisher_test_table, alternative = "greater") # 0.03823

ylimTop <- max(barplot_data)
tmp_barplot1 <- barplot(barplot_data, col = c("red", "blue"), 
                        beside = T, ylim = c(0, ylimTop + 5),
                        legend = c("CNV Breakpoints", "Random Breakpoints"), 
                        main = "Distance between TAD boundaries\nand CNV Breakpoints, RPMI-8226", 
                        axes=F, axisnames = F,
                        ylab = "Frequency of Distance",
                        xlab = "Distance", 
                        args.legend = list(x = "top", bty = 'n'))
segments(tmp_barplot1[2, ], test_random_tad_distance_1000_RPMI8226_counts_mean - 2*test_random_tad_distance_1000_RPMI8226_counts_sd, 
         tmp_barplot1[2, ], test_random_tad_distance_1000_RPMI8226_counts_mean + 2*test_random_tad_distance_1000_RPMI8226_counts_sd, 
         col = "black", lwd = 1.5)

arrows(tmp_barplot1[2, ], test_random_tad_distance_1000_RPMI8226_counts_mean - 2*test_random_tad_distance_1000_RPMI8226_counts_sd, 
       tmp_barplot1[2, ], test_random_tad_distance_1000_RPMI8226_counts_mean + 2*test_random_tad_distance_1000_RPMI8226_counts_sd, 
       lwd = 1.5, angle = 90,
       code = 3, length = 0.03)

axis( side = 1, at = tmp_barplot1[1, c(1, 6, 11, 16, 21, 26)], cex.axis = 0.7,
      label = c("0kb", "200kb", "400kb", "600kb", "800kb", "≥1000kb"), col = "black")
axis(side = 2)
dev.off()


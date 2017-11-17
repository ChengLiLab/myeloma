# A permutation experiment to test whether the cnv block boundaries is significantly overlapped with TAD boundaries 
# First select random sites, plot the frequency of overlapping sites numbers with TAD boundaries
# Then compare with the frequency of CNV sites (RPMI-8226 and U266 ) that overlapped with TAD boundraies
# 20160329
# wupz
# 1.3 read data of cnv
U266_cnv_40kb <- read.table( '/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_40kb/U266-merged.sorted.bam.dedup.bam_ratio.txt',
                                 header = T,
                                 stringsAsFactors = F)
# Change chromose X Y to 23, 24
U266_cnv_40kb$Chromosome[ U266_cnv_40kb$Chromosome == 'X'] <- 23
U266_cnv_40kb$Chromosome[ U266_cnv_40kb$Chromosome == 'Y'] <- 24
U266_cnv_block_40kb <- FindCNVBlock( freec_ratio = U266_cnv_40kb)
dim(U266_cnv_block_40kb)

select_random_sites_one_chrom <- function( all_sites, number_of_selected_sites) {
  select_id <- sample( x = 1:dim(all_sites)[1], size = number_of_selected_sites)
  return( all_sites[select_id, ])
}
# read insualtion score date
U266_is_40kb <- list()
for ( i in 1:23 ) {
  U266_is_40kb[[i]] <- read.table( paste( '/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/TAD_boundary/chr',
                                              i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation', sep = '' ),
                                       header = T,
                                       stringsAsFactors = F)
}
# read TAD data
U266_is_40kb_TAD <- list()
for ( i in 1:23 ) {
  U266_is_40kb_TAD[[i]] <- read.table( paste( '/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/TAD_boundary/chr',
                                                  i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation.boundaries', sep = '' ),
                                           header = T,
                                           stringsAsFactors = F)
}
length(U266_is_40kb_TAD)
sum(sapply(U266_is_40kb_TAD, dim)[1, ])

U266_cnv_block_40kb_tad_distance <- numeric()
for ( i in 1:23 ) {
  tmp_cnv_tad_distance <- nearestTADDistance( query_Sites = U266_cnv_block_40kb[U266_cnv_block_40kb[,1] == i, ],  subject_Sites =  U266_is_40kb_TAD[[i]])
  U266_cnv_block_40kb_tad_distance <- c(U266_cnv_block_40kb_tad_distance, tmp_cnv_tad_distance)
}
U266_cnv_block_40kb_numbers <- dim(U266_cnv_block_40kb)[1]
sum(U266_cnv_block_40kb_tad_distance == 0)/U266_cnv_block_40kb_numbers #  [1] 0.06976744        # 30/430
sum(U266_cnv_block_40kb_tad_distance <= 80000)/U266_cnv_block_40kb_numbers # [1] 0.2883721      # 124/430

# test for 1000 times
# test_random_tad_distance_1000: 40 kb
test_random_tad_distance_1000_U266 <- numeric()
for ( i in 1:1000) {
  test_random_tad_distance_1000_U266 <- cbind( test_random_tad_distance_1000_U266, 
                                          random_tad_distance( insulation_score = U266_is_40kb, TAD_insulation_score = U266_is_40kb_TAD, cnv_block_sites = U266_cnv_block_40kb) 
  )
}
png(filename = "../CNV_TAD_boundary/U266 Random regions overlapped with TAD Boundaries.png", width = 1024, height = 1024)
par( cex = 4, mgp = c(1.5, 0.5, 0), pch = '.', lwd = 4, mar = c(3, 4, 4, 3))
tmp_hist <- hist( apply(test_random_tad_distance_1000_U266, 2, function(x) sum(x<=80000)), breaks = 18, 
                  main = "Random regions overlapped\n with TAD Boundaries, U266", xlab = "Permutation times",
                  ylab = 'Frequency of Overlapped Regions', col = 'deepskyblue1', cex.axis = 0.8 )
arrows(x0 = tmp_hist$mids[11], y0 = 160, x1 = tmp_hist$mids[11], y1 = 0, col = 'red', lwd = 6) # x0 is the position of 124 which is the number of overlap regions
tmp_pval <- sum(apply(test_random_tad_distance_1000, 2, function(x) sum(x<=80000)) >= 124)/1000
text( x = tmp_hist$mids[11], y = 190, cex = 0.55, col = "red",
      labels = paste('CNV block boundaries \noverlapped with \nTAD Boundaries\nP value:', tmp_pval))
dev.off()

setwd("../CNV_TAD_boundary")
save.image('20160802.RData')
load("../CNV_TAD_boundary/20160329.RData")

png(filename = "../CNV_TAD_boundary/U266 Distance of CNV boundaries to TAD Boundaries.png", width = 1024, height = 1024)
par( cex = 3.5, mgp = c(2.5, 0.5, 0), pch = '.', lwd = 4, mar = c(5, 4, 4, 3))
tmp_coordinate <- barplot(table(U266_cnv_block_40kb_tad_distance), 
                          main = "Distance between TAD boundaries\n and CNV block boundaries, U266", 
                          ylab = "Frequency of Distance", col = "orangered",
                          xaxt = "n")
length(tmp_coordinate) # 56 positions
table(U266_cnv_block_40kb_tad_distance) # 56 distances
cbind(1:length(tmp_coordinate), tmp_coordinate, table(U266_cnv_block_40kb_tad_distance)) # found corresponding coordinate
axis(side = 1, at = tmp_coordinate[c(1, 13, 26, 37, 46, 51)], cex.axis = 0.65, 
     labels = paste( c(0, 0.5, 1, 1.5, 2, 20), "M" ) )

png(filename = "../CNV_TAD_boundary/U266 Distance of CNV boundaries to TAD Boundaries.png", width = 1024, height = 1024)
par( cex = 3, mgp = c(2, 0.5, 0), pch = '.', lwd = 4, mar = c(5, 4, 4, 3))
U266_cnv_block_40kb_tad_distance_counts <- table(cut(U266_cnv_block_40kb_tad_distance, breaks = c(0:25*40000, Inf)-1, include.lowest = T, right = F))
test_random_tad_distance_1000_U266_counts <- apply(X = test_random_tad_distance_1000_U266, 2, FUN = function(x) table(cut(x, breaks = c(0:25*40000, Inf)-1, include.lowest = T, right = F)))
test_random_tad_distance_1000_U266_counts_mean <- apply(test_random_tad_distance_1000_U266_counts, 1, mean)
test_random_tad_distance_1000_U266_counts_sd <- apply(test_random_tad_distance_1000_U266_counts, 1, sd)
barplot_data <- rbind(U266_cnv_block_40kb_tad_distance_counts, 
                      test_random_tad_distance_1000_U266_counts_mean)
# One Side Fisher Exact Test: 0.04292
fisher_test_table <- matrix(c(sum(U266_cnv_block_40kb_tad_distance_counts[1:3] ), 
                              sum(U266_cnv_block_40kb_tad_distance_counts[4:26] ), 
                              sum(test_random_tad_distance_1000_U266_counts_mean[1:3] ),
                              sum(test_random_tad_distance_1000_U266_counts_mean[4:26])), 
                            nrow = 2, byrow = T)
fisher.test(fisher_test_table, alternative = "greater")
ylimTop <- max(barplot_data)
tmp_barplot1 <- barplot(barplot_data, col = c("red", "blue"), 
                        beside = T, ylim = c(0, ylimTop + 5),
                        legend = c("CNV Breakpoints", "Romdom Sites"), 
                        main = "Distance between TAD boundaries\nand CNV Breakpoints, U266", 
                        axes=F, axisnames = F,
                        ylab = "Frequency of Distance",
                        xlab = "Distance", 
                        args.legend = list(x = "top", bty = 'n'))
segments(tmp_barplot1[2, ], test_random_tad_distance_1000_U266_counts_mean - 2*test_random_tad_distance_1000_U266_counts_sd, 
         tmp_barplot1[2, ], test_random_tad_distance_1000_U266_counts_mean + 2*test_random_tad_distance_1000_U266_counts_sd, 
         col = "black", lwd = 1.5)

arrows(tmp_barplot1[2, ], test_random_tad_distance_1000_U266_counts_mean - 2*test_random_tad_distance_1000_U266_counts_sd, 
       tmp_barplot1[2, ], test_random_tad_distance_1000_U266_counts_mean + 2*test_random_tad_distance_1000_U266_counts_sd, 
       lwd = 1.5, angle = 90,
       code = 3, length = 0.03)

axis( side = 1, at = tmp_barplot1[1, c(1, 6, 11, 16, 21, 26)], cex.axis = 0.7,
      label = c("0kb", "200kb", "400kb", "600kb", "800kb", "≥1000kb"), col = "black")
axis(side = 2)

dev.off()


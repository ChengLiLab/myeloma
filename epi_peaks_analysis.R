table(RPMI8226_cnv_40kb[,1])
cover.bed <- cbind(paste('chr', RPMI8226_cnv_40kb[,1], sep = ''), RPMI8226_cnv_40kb[,2], RPMI8226_cnv_40kb[,2] + 200000)
head(cover.bed)
write.table(cover.bed, file = '/lustre/user/liclab/wupz/dosageEffect/epigenome/data/cover.bed', sep = '\t', 
            append = F, quote = F, row.names = F, col.names = F)


peaks_1 <- read.table(file = '/lustre/user/liclab/jialm/DATA/publicDATA/chip/U266/H3K27ac/txt/U266_H3k27ac_summits.bed')
head(peaks_1)

near_cnv_peak_1 <- numeric()
near_cnv_peak_1_random <- numeric()
for ( i in 1:22 ) {
  tmp_vec <- nearestTADDistance(  peaks_1[peaks_1[, 1] == paste('chr', i, sep = ''), ], RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ])
  near_cnv_peak_1 <- c(near_cnv_peak_1, tmp_vec)
  tmp_random_number <- dim(RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ] )[1]
  tmp_random_site <- RPMI8226_is_40kb[[i]][ sample( x = 1:dim(RPMI8226_is_40kb[[i]])[1], size = tmp_random_number), ]
  tmp_nearest_distance <- nearestTADDistance( query_Sites = tmp_random_site, subject_Sites = RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ])
  near_cnv_peak_1_random <- c(near_cnv_peak_1_random, tmp_nearest_distance)
}
png(filename = 'peak_distance_h3k27ac.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1)
plot(density(near_cnv_peak_1),col = 'red', xlab = 'Distance', main = 'Distance of H3K27ac Peaks \nto CNV block Boundaries')
lines(density(near_cnv_peak_1_random) )
legend('topright', legend = c('CNV Boundaries', 'Random Sites'), title = 'RPMI-8226',
       col = c('red', 'black'), cex = 0.8, bty = 'n', lty = 1)
dev.off()
median(near_cnv_peak_1)
median(near_cnv_peak_1_random)

peaks_2 <- read.table(file = '/lustre/user/liclab/jialm/DATA/publicDATA/chip/U266/H3K27me3/txt/U266_H3k27me3_summits.bed')
near_cnv_peak_2 <- numeric()
near_cnv_peak_2_random <- numeric()
for ( i in 1:22 ) {
  tmp_vec <- nearestTADDistance(  peaks_2[peaks_2[, 1] == paste('chr', i, sep = ''), ], RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ])
  near_cnv_peak_2 <- c(near_cnv_peak_2, tmp_vec)
  tmp_random_number <- dim(RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ] )[1]
  tmp_random_site <- RPMI8226_is_40kb[[i]][ sample( x = 1:dim(RPMI8226_is_40kb[[i]])[1], size = tmp_random_number), ]
  tmp_nearest_distance <- nearestTADDistance( query_Sites = tmp_random_site, subject_Sites = RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ])
  near_cnv_peak_2_random <- c(near_cnv_peak_2_random, tmp_nearest_distance)
}
png(filename = 'peak_distance_H3K27me3.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1)
plot(density(near_cnv_peak_2),col = 'red', xlab = 'Distance', main = 'Distance of H3K27me3 Peaks \nto CNV block Boundaries')
lines(density(near_cnv_peak_2_random) )
legend('topright', legend = c('CNV Boundaries', 'Random Sites'), title = 'RPMI-8226',
       col = c('red', 'black'), cex = 0.8, bty = 'n', lty = 1)
dev.off()
median(near_cnv_peak_2)
median(near_cnv_peak_2_random)

peaks_3 <- read.table(file = '/lustre/user/liclab/jialm/DATA/publicDATA/chip/U266/H3K36me3/txt/U266_H3k36me3_summits.bed')
near_cnv_peak_3 <- numeric()
near_cnv_peak_3_random <- numeric()
for ( i in 1:22 ) {
  tmp_vec <- nearestTADDistance(  peaks_3[peaks_3[, 1] == paste('chr', i, sep = ''), ], RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ])
  near_cnv_peak_3 <- c(near_cnv_peak_3, tmp_vec)
  tmp_random_number <- dim(RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ] )[1]
  tmp_random_site <- RPMI8226_is_40kb[[i]][ sample( x = 1:dim(RPMI8226_is_40kb[[i]])[1], size = tmp_random_number), ]
  tmp_nearest_distance <- nearestTADDistance( query_Sites = tmp_random_site, subject_Sites = RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ])
  near_cnv_peak_3_random <- c(near_cnv_peak_3_random, tmp_nearest_distance)
}
png(filename = 'peak_distance_H3K36me3.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1)
plot(density(near_cnv_peak_3),col = 'red', xlab = 'Distance', main = 'Distance of H3K36me3 Peaks \nto CNV block Boundaries')
lines(density(near_cnv_peak_3_random) )
legend('topright', legend = c('CNV Boundaries', 'Random Sites'), title = 'RPMI-8226',
       col = c('red', 'black'), cex = 0.8, bty = 'n', lty = 1)
dev.off()
median(near_cnv_peak_3)
median(near_cnv_peak_3_random)
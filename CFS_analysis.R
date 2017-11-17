# CFS analysis
# 20160322
# wupz

# 1.0 read data
# 1.1 read CFS data
lym_cfs <- read.table("CFS_HG19.txt", skip = 3)
head(lym_cfs)
colnames(lym_cfs) <- c('chrom', 'start', 'end')

nearCFS <- numeric()
nearRandom <- numeric()
for ( i in 1:22 ) {
  tmp_vec <- nearestTADDistance(  RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ], lym_cfs[lym_cfs[, 1] == paste('chr',i,sep = ''), ])
  nearCFS <- c(nearCFS, tmp_vec)
  tmp_random_number <- dim(RPMI8226_cnv_block_40kb[RPMI8226_cnv_block_40kb[,1] == i, ] )[1]
  tmp_random_site <- RPMI8226_is_40kb[[i]][ sample( x = 1:dim(RPMI8226_is_40kb[[i]])[1], size = tmp_random_number), ]
  tmp_nearest_distance <- nearestTADDistance( query_Sites = tmp_random_site, subject_Sites =lym_cfs[lym_cfs[, 1] == paste('chr',i,sep = ''), ] )
  nearRandom <- c(nearRandom, tmp_nearest_distance)
}
png(filename = 'RPMI8226_CNV_boundaries_cfs_distance.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1)
plot(density(nearRandom),  xlab = 'Distance to CFSs', ylab = 'Probability Density', 
     main = 'Distance between CNV boundaris and CFSs\nCompared with random sites')
lines(density(nearCFS), col = 'red' )
legend('topright', legend = c('CNV Boundaries', 'Random Sites'), title = 'RPMI-8226',
       col = c('red', 'black'), cex = 0.8, bty = 'n', lty = 1)
dev.off()
sum( nearRandom == 0 )
sum( nearRandom <= 500000 )
median(nearRandom )
length(nearRandom)
sum( nearCFS == 0 )
sum( nearCFS <= 500000 )
median(nearCFS )

one_cfs_test <- function(cnv_block ) {
  nearRandom <- numeric()
  for ( i in 1:22 ) {
    tmp_random_number <- dim(cnv_block[cnv_block[,1] == i, ] )[1]
    tmp_random_site <- RPMI8226_is_40kb[[i]][ sample( x = 1:dim(RPMI8226_is_40kb[[i]])[1], size = tmp_random_number), ]
    tmp_nearest_distance <- nearestTADDistance( query_Sites = tmp_random_site, subject_Sites = lym_cfs[lym_cfs[, 1] == paste('chr',i,sep = ''), ] )
    nearRandom <- c(nearRandom, tmp_nearest_distance)
  }
  nearRandom_zeros <- sum(nearRandom == 0 )
  nearRandom_500kb <- sum(nearRandom <= 500000 )
  nearRandom_median <- median( nearRandom)
  return(c(nearRandom_zeros, nearRandom_500kb, nearRandom_median))
}

RPMI8226_random_test <- numeric()
for ( i in 1:1000) {
  tmp_random_test <- one_cfs_test( cnv_block = RPMI8226_cnv_block_40kb)
  RPMI8226_random_test <- rbind(RPMI8226_random_test, tmp_random_test)
}
png(filename = 'RPMI8226_CNV_boudaries_within_CFS.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1)
tmp_p <- sum(RPMI8226_random_test[, 1] <= sum(nearCFS == 0))/1000
plot(density(RPMI8226_random_test[, 1]), xlab =  'Proportion',
     main = paste('Proportion of sites within CFSs\n p value: ', tmp_p),'')
abline( v = sum(nearCFS == 0 ), col = 'red')
legend('topright', col = c('red', 'black'), legend = c('CNV Boundaries', 'Random Sites'), bty = 'n',
       title = 'RPMI-8226', lty = 1)
dev.off()

png(filename = 'RPMI8226_CNV_boudaries_within_500kb_CFS.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1)
tmp_p <- sum(RPMI8226_random_test[, 2] <= sum(nearCFS <= 500000))/1000
plot(density(RPMI8226_random_test[, 2]), xlab =  'Proportion',
     main = paste('Proportion of sites within 500kb around CFSs\n p value: ', tmp_p),'')
abline( v = sum(nearCFS <= 500000), col = 'red')
legend('topright', col = c('red', 'black'), legend = c('CNV Boundaries', 'Random Sites'), bty = 'n',
       title = 'RPMI-8226', lty = 1)
dev.off()

png(filename = 'RPMI8226_distance_of_CNV_boudaries_and_CFS.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1)
tmp_p <- sum(RPMI8226_random_test[, 3] >= median(nearCFS) )/1000
plot(density(RPMI8226_random_test[, 3]), xlab =  'Distance',
     main = paste('Median distance of sites with CFSs\n p value: ', tmp_p),'')
abline( v = median(nearCFS), col = 'red')
legend('topright', col = c('red', 'black'), legend = c('CNV Boundaries', 'Random Sites'), bty = 'n',
       title = 'RPMI-8226', lty = 1)
dev.off()


# 2. U266
U266_nearCFS <- numeric()
U266_nearRandom <- numeric()
for ( i in 1:22 ) {
  tmp_vec <- nearestTADDistance(  U266_cnv_block_40kb[U266_cnv_block_40kb[,1] == i, ], lym_cfs[lym_cfs[, 1] == paste('chr',i,sep = ''), ])
  U266_nearCFS <- c(U266_nearCFS, tmp_vec)
  tmp_random_number <- dim(U266_cnv_block_40kb[U266_cnv_block_40kb[,1] == i, ] )[1]
  tmp_random_site <- U266_is_40kb[[i]][ sample( x = 1:dim(U266_is_40kb[[i]])[1], size = tmp_random_number), ]
  tmp_nearest_distance <- nearestTADDistance( query_Sites = tmp_random_site, subject_Sites =lym_cfs[lym_cfs[, 1] == paste('chr',i,sep = ''), ] )
  U266_nearRandom <- c(U266_nearRandom, tmp_nearest_distance)
}
png(filename = 'U266_CNV_boundaries_cfs_distance.png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1)
plot(density(U266_nearRandom), col = 'black',  xlab = 'Distance to CFSs', ylab = 'Probability Density', 
     main = 'Distance between CNV boundaris and CFSs\nCompared with random sites' ) 
lines(density(U266_nearCFS), col = 'red' )
legend('topright', legend = c('CNV Boundaries', 'Random Sites'), title = 'U266',
       col = c('red', 'black'), cex = 0.8, bty = 'n', lty = 1)
dev.off()

sum( U266_nearRandom == 0 )
sum( U266_nearRandom <= 500000 )
median(U266_nearRandom )
length(nearRandom )
sum( U266_nearCFS == 0 )
sum( U266_nearCFS <= 500000 )
median(U266_nearCFS )

U266_random_test <- numeric()
for ( i in 1:1000) {
  tmp_random_test <- one_cfs_test( cnv_block = U266_cnv_block_40kb)
  U266_random_test <- rbind(U266_random_test, tmp_random_test)
}

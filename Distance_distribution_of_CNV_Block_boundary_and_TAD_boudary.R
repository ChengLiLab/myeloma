# Distance distribution of CNV Block boundary and TAD boudary
# 20160113
# wupz

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

# 2.2 get the distance of each CNV block boudnary
drawDistance <- function(cnv_sites, subject_sites_list, all_sites, ...) {
  CNV_block_nearest_distance <- numeric()
  tamdom_nearest_distance <- numeric()
  for ( i in 1:23) {
    tmp_nearest_distance <- nearestTADDistance( query_Sites = cnv_sites[cnv_sites[,1] == i, ], 
                                                subject_Sites =  subject_sites_list[[i]])
    CNV_block_nearest_distance <- c(CNV_block_nearest_distance, tmp_nearest_distance)
    tmp_random_number <- dim(cnv_sites[cnv_sites[,1] == i, ])[1]  # for each chromosome random select some sites with the same number of TADs
    tmp_random_site <- all_sites[[i]][ sample( x = 1:dim(all_sites[[i]])[1], size = tmp_random_number), ]
    tmp_nearest_distance <- nearestTADDistance( query_Sites = tmp_random_site, subject_Sites = subject_sites_list[[i]])
    tamdom_nearest_distance <- c(tamdom_nearest_distance, tmp_nearest_distance)
  }
  density_tamdom <- density(tamdom_nearest_distance )
  density_query <- density(CNV_block_nearest_distance )
  ylim = c(0, max(density_tamdom$y, density_query$y)) 
  plot(density(tamdom_nearest_distance ), col = 'black', ...)
  lines(density(CNV_block_nearest_distance), col = 'red' )
}

png(filename = 'Distance Distribution of sites from TAD boudary RPMI8226.png', height = 1024, width = 1024)
par(cex = 1.5, mgp = c(1, 0.2, 0), lwd = 3, tcl =-0.1)
drawDistance( RPMI8226_cnv_block_40kb, subject_sites_list = RPMI8226_is_40kb_TAD, RPMI8226_is_40kb, xlim = c(0, 500000),
              main = 'Distance Distribution of sites from TAD boudary, RPMI8226')
legend( 'topright', legend = c('CNV block boudary', 'Random'), col = c('Red', 'black'), lty = 1)
dev.off()

png(filename = 'Distance Distribution of sites from TAD boudary U266', height = 1024, width = 1024)
par(cex = 1.5, mgp = c(1, 0.2, 0), lwd = 3, tcl =-0.1)
drawDistance( U266_cnv_block_40kb, subject_sites_list = U266_is_40kb_TAD, U266_is_40kb, xlim = c(0, 500000),
              main = 'Distance Distribution of sites from TAD boudary, U266')
legend( 'topright', legend = c('CNV block boudary', 'Random'), col = c('Red', 'black'), lty = 1)
dev.off()

png(filename = 'Distance Distribution of sites from TAD boudary GM12878', height = 1024, width = 1024)
par(cex = 1.5, mgp = c(1, 0.2, 0), lwd = 3, tcl =-0.1)
drawDistance( U266_cnv_block_40kb, subject_sites_list = GM12878_is_10kb_TAD, GM12878_is_10kb, xlim = c(0, 500000),
              main = 'Distance Distribution of sites from TAD boudary, GM12878')
legend( 'topright', legend = c('CNV block boudary', 'Random'), col = c('Red', 'black'), lty = 1)
dev.off()



CNV_block_nearest_distance <- numeric()
tamdom_nearest_distance <- numeric()
for ( i in 1:23) {
  tmp_nearest_distance <- nearestTADDistance( query_Sites = U266_cnv_block_40kb[U266_cnv_block_40kb[,1] == i, ], 
                                              subject_Sites =  GM12878_is_10kb_TAD[[i]])
  CNV_block_nearest_distance <- c(CNV_block_nearest_distance, tmp_nearest_distance)
  tmp_random_number <- dim(U266_cnv_block_40kb[U266_cnv_block_40kb[,1] == i, ])[1]  # for each chromosome random select some sites with the same number of TADs
  tmp_random_site <- GM12878_is_10kb[[i]][ sample( x = 1:dim(U266_is_40kb[[i]])[1], size = tmp_random_number), ]
  tmp_nearest_distance <- nearestTADDistance( query_Sites = tmp_random_site, subject_Sites = GM12878_is_10kb_TAD[[i]])
  tamdom_nearest_distance <- c(tamdom_nearest_distance, tmp_nearest_distance)
}

plot(density(tamdom_nearest_distance ), xlim = c(0, 1000000))
lines(density(CNV_block_nearest_distance), col = 'red')
median(tamdom_nearest_distance)
median(CNV_block_nearest_distance)
t.test(tamdom_nearest_distance, CNV_block_nearest_distance)

nearTAD_8226_U266 <- numeric()
for ( i in 1:23 ) {
  tmp_nearest_distance <- nearestTADDistance( query_Sites = RPMI8226_is_40kb_TAD[[i]], 
                                              subject_Sites =  U266_is_40kb_TAD[[i]])
  nearTAD_8226_U266 <- c(nearTAD_8226_U266, tmp_nearest_distance)
}
png(filename = 'The Shifts of TAD Boundaries (kb).png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1)
plot(table(nearTAD_8226_U266), xlab = 'The Shifts of TAD Boundaries ( bp)', ylab = 'RPMI-8226 v.s. U266', type = 'l',
     col = 'red', axes = FALSE, main = 'The Shifts of TAD Boundaries\n between RPMI-8226 and U266')
axis(2)
axis(1)
dev.off()
write.table( x = table(nearTAD_8226_U266), file = 'The Shifts of TAD Boundaries.txt')
getwd()

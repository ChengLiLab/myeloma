# analysis cnv boudary and TAD boundary by using public data
# 20160302
# wupz

# source some scripts
source("../scripts/Sites_IS.R")
source("../scripts/LocalIS.R")

# 1. read data
# 1.1 Database of Genomic Variants
CNV_dgv <- read.table( file = '/lustre/user/liclab/wupz/dosageEffect/rawData/GRCh37_hg19_variants_2014-10-16.txt',sep = '\t', 
                       header = T, stringsAsFactors = F)
head(CNV_dgv)
dim(CNV_dgv)

# 1.2 filter some data
# 1.2.1 keep large cnv (>= 500kb)
filter_id_1 <- ( CNV_dgv[, 'end'] - CNV_dgv[, 'start'] ) >= 500000
# 1.2.2 keep sample rate >= 0.1
filter_id_2 <- ( CNV_dgv[, 'observedgains'] - CNV_dgv[, 'observedlosses'] )/CNV_dgv[, 'samplesize'] >= 0.1
filter_id_2 <- filter_id_2 & CNV_dgv[, 'chr'] != "Y"
sum(filter_id_1)
sum(filter_id_1 & filter_id_2)
head(CNV_dgv[filter_id_1 & filter_id_2, ])
table(CNV_dgv[filter_id_1 & filter_id_2, 'chr'])
CNV_dgv[filter_id_1 & filter_id_2, 'chr'][CNV_dgv[filter_id_1 & filter_id_2, 'chr'] == 'X' ] <- 23 
CNV_dgv[filter_id_1 & filter_id_2, 'chr'][CNV_dgv[filter_id_1 & filter_id_2, 'chr'] == 'Y' ] <- 24
# GM12878_is_40kb
GM12878_is_40kb <- list()
for ( i in 1:23 ) {
  GM12878_is_40kb[[i]] <- read.table( paste('/lustre/user/liclab/lirf/Project/hic/2014/gm12878_hind3_dCTP_296mb/resolution_40k/cis/TAD_boundary/chr',
                                            i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation', sep = ''),
                                      header = T,
                                      stringsAsFactors = F)
}

# 1.2 read data of TAD boundary
GM12878_is_40kb_TAD <- list()
for ( i in  c(1:23)) {
  GM12878_is_40kb_TAD[[i]] <- read.table( paste('/lustre/user/liclab/lirf/Project/hic/2014/gm12878_hind3_dCTP_296mb/resolution_40k/cis/TAD_boundary/chr',
                                                i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation.boundaries', sep = ''),
                                          header = T,
                                          stringsAsFactors = F)
}
dim(GM12878_is_40kb_TAD[[3]])

# 2.1 for each cnv from dgv, find the is around the boundary, and the same with tandom sites， 10kb 
CNV_dgv_InsulationScore_GM12878  <- numeric()
random_insulationScore_GM12878 <- numeric()
for ( i in names(table(CNV_dgv[filter_id_1 , 'chr'])) ) {
  tmp_sites <-  CNV_dgv[filter_id_1 , 2:4]
  tmp_cnv_block_ins <- Sites_IS( Sites = tmp_sites[tmp_sites[,1] == i, ], insulation_score = GM12878_is_40kb[[i]], flank = 100)
  CNV_dgv_InsulationScore_GM12878 <- rbind(CNV_dgv_InsulationScore_GM12878, tmp_cnv_block_ins)
  # tandom sites 
  tmp_random_number <- dim(tmp_sites)[1] # for each chromosome random select some sites with the same number of TADs
  tmp_random_site <- GM12878_is_40kb[[i]][ sample( x = 1:dim(GM12878_is_40kb[[i]])[1], size = tmp_random_number), ]
  tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score = GM12878_is_40kb[[i]], flank = 100)
  random_insulationScore_GM12878   <- rbind(random_insulationScore_GM12878  ,  tmp_random_ins)
}
median_ins_CNV_dgv <- apply(CNV_dgv_InsulationScore_GM12878, 2, function(x) median( as.numeric(x), na.rm = T) )
median_random_insulationScore_GM12878 <- apply(random_insulationScore_GM12878, 2, function(x) median( as.numeric(x), na.rm = T) )

# 2.2 for each cnv from dgv, find the is around the boundary, and the same with tandom sites， 40kb 
CNV_dgv_InsulationScore_GM12878  <- numeric()
random_insulationScore_G M12878 <- numeric()
TAD_InsulationScore_GM12878_40kb <- numeric()
tmp_sites <-  as.matrix( CNV_dgv[filter_id_1 & filter_id_2, 2:4] )
for ( i in names(table(CNV_dgv[filter_id_1 & filter_id_2, 'chr'])) ) {
  if (table(CNV_dgv[filter_id_1 & filter_id_2, 'chr'])[i] > 1 ) {
    tmp_cnv_block_ins <- Sites_IS( Sites = tmp_sites[tmp_sites[,1] == i, ], insulation_score = GM12878_is_40kb[[ as.numeric(i)]], flank = 25)
    CNV_dgv_InsulationScore_GM12878 <- rbind(CNV_dgv_InsulationScore_GM12878, tmp_cnv_block_ins)
    # tandom sites 
    tmp_random_number <- dim(tmp_sites)[1] # for each chromosome random select some sites with the same number of TADs
    tmp_random_site <- GM12878_is_40kb[[as.numeric(i)]][ sample( x = 1:dim(GM12878_is_40kb[[as.numeric(i)]])[1], size = tmp_random_number), ]
    tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score = GM12878_is_40kb[[as.numeric(i)]], flank = 25)
    random_insulationScore_GM12878   <- rbind(random_insulationScore_GM12878  ,  tmp_random_ins)
    # GM12878 TAD boundary
    print(dim(GM12878_is_40kb_TAD[[ as.numeric(i)]]))
    tmp_TAD_ins <- Sites_IS( Sites = GM12878_is_40kb_TAD[[ as.numeric(i)]], insulation_score = GM12878_is_40kb[[ as.numeric(i)]], flank = 25) 
    TAD_InsulationScore_GM12878_40kb <- rbind(TAD_InsulationScore_GM12878_40kb, tmp_TAD_ins)
  }
}
median_ins_CNV_dgv <- apply(CNV_dgv_InsulationScore_GM12878, 2, function(x) median( as.numeric(x), na.rm = T) )
median_random_insulationScore_GM12878 <- apply(random_insulationScore_GM12878, 2, function(x) median( as.numeric(x), na.rm = T) )
GM12878_median_TAD_ins_40kb <- apply(TAD_InsulationScore_GM12878_40kb, 2, function(x) median( as.numeric(x), na.rm = T) )


# 2.2 draw the insulation score around the cnv sites
png(filename = 'Human genome cnv boundary and TAD boundary (GM12878).png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1, 0.2, 0), lwd = 3, tcl =-0.1)
plot( GM12878_median_TAD_ins_40kb, type = 'l', 
      xlab = 'Genomic Regions', 
      ylab = 'Median Insulation Score',
      main = 'Association between CNV Breakpoints\n and TAD boundaris',
      ylim = c(min(median_ins_CNV_dgv, GM12878_median_TAD_ins_40kb,median_random_insulationScore_GM12878, na.rm = T), 
               max(median_ins_CNV_dgv, GM12878_median_TAD_ins_40kb,median_random_insulationScore_GM12878, na.rm = T) ),
      col = "black",
      axes = FALSE)
lines( median_random_insulationScore_GM12878, col = "blue" )
lines(  median_ins_CNV_dgv, col = "red")

axis(side = 1, at = c(1, 14, 26, 38, 51), labels = c('-1Mb', '-480kb','0kb','480kb', '1Mb' ))
axis(side = 2)
legend('bottomright', bty = 'n', cex = 0.65, pch = 19,
       legend = c('Human Genome\nCNVs(>=500kb)', 'TAD Boudaries', 'Random Breakpoints'),
       col = c('red', 'black', 'blue' ) )
dev.off()


# 2.3 for each cnv from cnvd, find the is around the boundary, and the same with tandom sites 
CNV_CNVD_InsulationScore_GM12878  <- numeric()
random_CNVD_insulationScore_GM12878 <- numeric()
for ( i in names(table(CNV_CNVD[filter_id_3 & filter_id_4 & filter_id_5 , 'chr'])) ) {
  tmp_sites <-  CNV_CNVD[filter_id_3 & filter_id_4 & filter_id_5 , 3:5]
  tmp_sites <- tmp_sites[tmp_sites[,1] == i, ]
  tmp_cnv_block_ins <- Sites_IS( Sites = tmp_sites, insulation_score = GM12878_is_40kb[[i]], flank = 100)
  CNV_CNVD_InsulationScore_GM12878 <- rbind(CNV_CNVD_InsulationScore_GM12878, tmp_cnv_block_ins)
  # tandom sites 
  tmp_random_number <- dim(tmp_sites)[1] # for each chromosome random select some sites with the same number of TADs
  tmp_random_site <- GM12878_is_40kb[[i]][ sample( x = 1:dim(GM12878_is_40kb[[i]])[1], size = tmp_random_number), ]
  tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score = GM12878_is_40kb[[i]], flank = 100)
  random_CNVD_insulationScore_GM12878   <- rbind(random_CNVD_insulationScore_GM12878  ,  tmp_random_ins)
}
median_ins_CNV_CNVD <- apply(CNV_CNVD_InsulationScore_GM12878, 2, function(x) median( as.numeric(x), na.rm = T) )
median_random_insulationScore_GM12878_CNVD <- apply(random_CNVD_insulationScore_GM12878, 2, function(x) median( as.numeric(x), na.rm = T) )

# 2.4 draw the insulation score around the cnv sites
png(filename = 'Human disease (Multiple myeloma) cnv boundary and TAD boundary (GM12878).png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1)
plot( GM12878_median_TAD_ins, type = 'l', 
      xlab = 'Genomic Sites', 
      ylab = 'Median Insulation Score',
      main = 'CNV Block Boundaries and TAD boundaris',
      ylim = c(min(median_ins_CNV_CNVD, GM12878_median_TAD_ins,median_random_insulationScore_GM12878_CNVD, na.rm = T), 
               max(median_ins_CNV_CNVD, GM12878_median_TAD_ins,median_random_insulationScore_GM12878_CNVD, na.rm = T) ),
      col = "black",
      axes = FALSE)
lines( median_random_insulationScore_GM12878_CNVD, col = "grey" )
lines(  median_ins_CNV_CNVD, col = "red")
axis(side = 1, at = c(1, 50, 101, 151, 201), labels = c('-1Mb', '-500kb','0kb','500kb', '1Mb' ))
axis(side = 2)
legend('bottomright', bty = 'n', cex = 0.65,
       legend = c('TAD Boudaries', 'Random Sites', 'CNV block boundaries\n(Multiple Myeloma)'),
       col = c('black', 'grey', 'red'), lty = 1)
dev.off()


tmp_plotBoundaries <- function( tmp_vec ) {
  # 2.3 for each cnv from cnvd, find the is around the boundary, and the same with tandom sites
  filter_id_5 <- CNV_CNVD[, 'disease'] == tmp_vec
  CNV_CNVD_InsulationScore_GM12878  <- numeric()
  random_CNVD_insulationScore_GM12878 <- numeric()
  for ( i in names(table(CNV_CNVD[filter_id_3 & filter_id_4 & filter_id_5 , 'chr'])) ) {
    tmp_sites <-  CNV_CNVD[filter_id_3 & filter_id_4 & filter_id_5 , 3:5]
    tmp_sites <- tmp_sites[tmp_sites[,1] == i, ]
    tmp_cnv_block_ins <- Sites_IS( Sites = tmp_sites, insulation_score = GM12878_is_40kb[[i]], flank = 100)
    CNV_CNVD_InsulationScore_GM12878 <- rbind(CNV_CNVD_InsulationScore_GM12878, tmp_cnv_block_ins)
    # tandom sites 
    tmp_random_number <- dim(tmp_sites)[1] # for each chromosome random select some sites with the same number of TADs
    tmp_random_site <- GM12878_is_40kb[[i]][ sample( x = 1:dim(GM12878_is_40kb[[i]])[1], size = tmp_random_number), ]
    tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score = GM12878_is_40kb[[i]], flank = 100)
    random_CNVD_insulationScore_GM12878   <- rbind(random_CNVD_insulationScore_GM12878  ,  tmp_random_ins)
  }
  median_ins_CNV_CNVD <- apply(CNV_CNVD_InsulationScore_GM12878, 2, function(x) median( as.numeric(x), na.rm = T) )
  median_random_insulationScore_GM12878_CNVD <- apply(random_CNVD_insulationScore_GM12878, 2, function(x) median( as.numeric(x), na.rm = T) )
  
  # 2.4 draw the insulation score around the cnv sites
  png(filename = paste('Human disease (', tmp_vec, ') cnv boundary and TAD boundary (GM12878).png', sep = ''), width = 1024, height = 1024)
  par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1)
  plot( GM12878_median_TAD_ins, type = 'l', 
        xlab = 'Genomic Sites', 
        ylab = 'Median Insulation Score',
        main = 'CNV Block Boundaries and TAD boundaris',
        ylim = c(min(median_ins_CNV_CNVD, GM12878_median_TAD_ins,median_random_insulationScore_GM12878_CNVD, na.rm = T), 
                 max(median_ins_CNV_CNVD, GM12878_median_TAD_ins,median_random_insulationScore_GM12878_CNVD, na.rm = T) ),
        col = "black",
        axes = FALSE)
  lines( median_random_insulationScore_GM12878_CNVD, col = "grey" )
  lines(  median_ins_CNV_CNVD, col = "red")
  axis(side = 1, at = c(1, 50, 101, 151, 201), labels = c('-1Mb', '-500kb','0kb','500kb', '1Mb' ))
  axis(side = 2)
  legend('bottomright', bty = 'n', cex = 0.65,
         legend = c('TAD Boudaries', 'Random Sites', paste('CNV block boundaries\n(', tmp_vec, ')', sep = '') ),
         col = c('black', 'grey', 'red'), lty = 1)
  dev.off()
}

tmp_plotBoundaries('Acute myeloid leukemia')
tmp_plotBoundaries('Burkitts lymphoma')
tmp_plotBoundaries('Basal cell lymphoma')
tmp_plotBoundaries('Diffuse large b-cell lymphoma')
tmp_plotBoundaries('Chronic lymphocytic leukemia')
tmp_plotBoundaries('Hodgkins lymphoma')
tmp_plotBoundaries('small cell lung cancer')
tmp_plotBoundaries('Colorectal cancer')

for ( i in names(table(CNV_CNVD[filter_id_3 & filter_id_4, 'disease'] )[table(CNV_CNVD[filter_id_3 & filter_id_4, 'disease'] ) >= 100 ]) ) {
  print(i)
  try (tmp_plotBoundaries( tmp_vec = i))
}

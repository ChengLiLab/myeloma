# analysis cnv boudary and TAD boundary by using public data
# 20160302
# wupz

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
sum(filter_id_1)
sum(filter_id_1 & filter_id_2)
head(CNV_dgv[filter_id_1 & filter_id_2, ])
table(CNV_dgv[filter_id_1 & filter_id_2, 'chr'])
CNV_dgv[filter_id_1 & filter_id_2, 'chr'][CNV_dgv[filter_id_1 & filter_id_2, 'chr'] == 'X' ] <- 23 

# 1.3 CNVD database: disease CNV
CNV_CNVD <- read.table( file = '/lustre/user/liclab/wupz/dosageEffect/rawData/CNVD_All_Result_20160303.txt',sep = '\t', 
                        stringsAsFactors = F)
colnames(CNV_CNVD ) <- c( 'id', 'species', 'chr', 'start', 'end', 'band', 'type', 'genes', 'disease', 'methods', 'sample_size', 'rate', 'pubmed_id')
dim(CNV_CNVD)
head(CNV_CNVD)
# 1.4 filtering some data 
filter_id_3 <- (CNV_CNVD[, 'end'] - CNV_CNVD[, 'start'] ) >= 500000
sum(filter_id_3)
filter_id_4 <- CNV_CNVD[, 'species'] == 'Homo sapiens'
sum(filter_id_4)
table(CNV_CNVD[filter_id_3 & filter_id_4, 'disease'] )[table(CNV_CNVD[filter_id_3 & filter_id_4, 'disease'] ) >= 100]
filter_id_5 <- CNV_CNVD[, 'disease'] == 'Multiple myeloma'
# CNV_CNVD[, 'disease'] == 'Basal cell lymphoma' | CNV_CNVD[, 'disease'] == 'Acute myeloid leukemia' | CNV_CNVD[, 'disease'] == 'Burkitts lymphoma' | 
sum(filter_id_5)
sum(filter_id_3 & filter_id_4 & filter_id_5)

head( CNV_CNVD[filter_id_3 & filter_id_4, ])

# 2.1 for each cnv from dgv, find the is around the boundary, and the same with tandom sites 
CNV_dgv_InsulationScore_GM12878  <- numeric()
random_insulationScore_GM12878 <- numeric()
for ( i in names(table(CNV_dgv[filter_id_1 , 'chr'])) ) {
  tmp_sites <-  CNV_dgv[filter_id_1 , 2:4]
  tmp_cnv_block_ins <- Sites_IS( Sites = tmp_sites[tmp_sites[,1] == i, ], insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
  CNV_dgv_InsulationScore_GM12878 <- rbind(CNV_dgv_InsulationScore_GM12878, tmp_cnv_block_ins)
  # tandom sites 
  tmp_random_number <- dim(tmp_sites)[1] # for each chromosome random select some sites with the same number of TADs
  tmp_random_site <- GM12878_is_10kb[[i]][ sample( x = 1:dim(GM12878_is_10kb[[i]])[1], size = tmp_random_number), ]
  tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
  random_insulationScore_GM12878   <- rbind(random_insulationScore_GM12878  ,  tmp_random_ins)
}
median_ins_CNV_dgv <- apply(CNV_dgv_InsulationScore_GM12878, 2, function(x) median( as.numeric(x), na.rm = T) )
median_random_insulationScore_GM12878 <- apply(random_insulationScore_GM12878, 2, function(x) median( as.numeric(x), na.rm = T) )

# 2.2 draw the insulation score around the cnv sites
png(filename = 'Human genome cnv boundary and TAD boundary (GM12878).png', width = 1024, height = 1024)
plot( GM12878_median_TAD_ins, type = 'l', 
      xlab = 'Genomic Sites', 
      ylab = 'Median Insulation Score',
      main = 'CNV Block Boundaries and TAD boundaris',
      ylim = c(min(median_ins_CNV_dgv, GM12878_median_TAD_ins,median_random_insulationScore_GM12878, na.rm = T), 
               max(median_ins_CNV_dgv, GM12878_median_TAD_ins,median_random_insulationScore_GM12878, na.rm = T) ),
      col = "black",
      axes = FALSE)
lines( median_random_insulationScore_GM12878, col = "grey" )
lines(  median_ins_CNV_dgv, col = "red")

axis(side = 1, at = c(1, 50, 101, 151, 201), labels = c('-1Mb', '-500kb','0kb','500kb', '1Mb' ))

axis(side = 2)
legend('bottomright', bty = 'n', cex = 0.65,
       legend = c('TAD Boudaries', 'Random Sites', 'Human Genome CNVs(>=500kb)'),
       col = c('black', 'grey', 'red'), lty = 1)
dev.off()


# 2.3 for each cnv from cnvd, find the is around the boundary, and the same with tandom sites 
CNV_CNVD_InsulationScore_GM12878  <- numeric()
random_CNVD_insulationScore_GM12878 <- numeric()
for ( i in names(table(CNV_CNVD[filter_id_3 & filter_id_4 & filter_id_5 , 'chr'])) ) {
  tmp_sites <-  CNV_CNVD[filter_id_3 & filter_id_4 & filter_id_5 , 3:5]
  tmp_sites <- tmp_sites[tmp_sites[,1] == i, ]
  tmp_cnv_block_ins <- Sites_IS( Sites = tmp_sites, insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
  CNV_CNVD_InsulationScore_GM12878 <- rbind(CNV_CNVD_InsulationScore_GM12878, tmp_cnv_block_ins)
  # tandom sites 
  tmp_random_number <- dim(tmp_sites)[1] # for each chromosome random select some sites with the same number of TADs
  tmp_random_site <- GM12878_is_10kb[[i]][ sample( x = 1:dim(GM12878_is_10kb[[i]])[1], size = tmp_random_number), ]
  tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
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
    tmp_cnv_block_ins <- Sites_IS( Sites = tmp_sites, insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
    CNV_CNVD_InsulationScore_GM12878 <- rbind(CNV_CNVD_InsulationScore_GM12878, tmp_cnv_block_ins)
    # tandom sites 
    tmp_random_number <- dim(tmp_sites)[1] # for each chromosome random select some sites with the same number of TADs
    tmp_random_site <- GM12878_is_10kb[[i]][ sample( x = 1:dim(GM12878_is_10kb[[i]])[1], size = tmp_random_number), ]
    tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
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

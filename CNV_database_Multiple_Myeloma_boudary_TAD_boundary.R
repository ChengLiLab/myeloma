# analysis cnv boudary and TAD boundary by using public data
# 20160802
# wupz

# 1.1 CNVD database: disease CNV
CNV_CNVD <- read.table( file = '/lustre/user/liclab/wupz/dosageEffect/rawData/CNVD_HS_MM.txt',sep = '\t', 
                        stringsAsFactors = F)
colnames(CNV_CNVD ) <- c( 'id', 'species', 'chr', 'start', 'end', 'band', 'type', 'genes', 'disease', 'methods', 'sample_size', 'rate', 'pubmed_id')
dim(CNV_CNVD) # [1] 302  13
head(CNV_CNVD)
# 1.2 filtering some data 
filter_id_3 <- (CNV_CNVD[, 'end'] - CNV_CNVD[, 'start'] ) >= 500000
sum(filter_id_3) # [1] 277
filter_id_4 <- CNV_CNVD$rate >= 0.1
sum(filter_id_4) # [1] 183
#table(CNV_CNVD[filter_id_3 & filter_id_4, 'disease' ] )[table(CNV_CNVD[filter_id_3 & filter_id_4, 'disease'] ) >= 100]
#filter_id_5 <- CNV_CNVD[, 'disease'] == 'Multiple myeloma'
# CNV_CNVD[, 'disease'] == 'Basal cell lymphoma' | CNV_CNVD[, 'disease'] == 'Acute myeloid leukemia' | CNV_CNVD[, 'disease'] == 'Burkitts lymphoma' | 
#sum(filter_id_5) # [1] 302
sum(filter_id_3 & filter_id_4)
CNV_CNVD_MM_filter <- CNV_CNVD[filter_id_3 & filter_id_4, ]
CNV_CNVD_MM_filter_breakpoints <- 

# 2.1 for each cnv from cnvd, find the is around the boundary, and the same with tandom sites 
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
png(filename = '../CNV_TAD_boundary/Human disease (Multiple myeloma) cnv boundary and TAD boundary (GM12878).png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1, mar = c(3,3,3,2))
plot( median_ins_CNV_CNVD, type = 'l', 
      xlab = 'Genomic Regions', 
      ylab = 'Median Insulation Score',
      main = 'CNV Breakpoints in Cancers \nand TAD boundaris in GM12878',
      ylim = c(min(median_ins_CNV_CNVD, GM12878_median_TAD_ins,median_random_insulationScore_GM12878_CNVD, na.rm = T), 
               max(median_ins_CNV_CNVD, GM12878_median_TAD_ins,median_random_insulationScore_GM12878_CNVD, na.rm = T) ),
      col = "red",
      lwd = 5, 
      xlim = c(0, 250), 
      axes = FALSE)
lines( GM12878_median_TAD_ins, lwd = 5, col = "black" )
lines(median_random_insulationScore_GM12878_CNVD, lwd = 5, col = "blue")
axis(side = 1, at = c(1, 50, 101, 151, 201), labels = c('-1000kb', '-500kb','0kb','500kb', '1000kb' ))
axis(side = 2)
legend('bottomright', bty = 'n', cex = 0.75, 
       legend = c('CNV Breakpoints\n(Multiple Myeloma)', 'TAD Boudaries', 'Random Sites'),
       col = c( 'red', 'black', 'blue'), pch = 19)
dev.off()


# select a disease : 
filter_id_6 <- CNV_CNVD[, 'disease'] == 'High-grade prostatic intraepithelial neoplasia'
# CNV_CNVD[, 'disease'] == 'Basal cell lymphoma' | CNV_CNVD[, 'disease'] == 'Acute myeloid leukemia' | CNV_CNVD[, 'disease'] == 'Burkitts lymphoma' | 
sum(filter_id_6) # [1] 208
sum(filter_id_3 & filter_id_4 & filter_id_6)
head( CNV_CNVD[filter_id_3 & filter_id_4, ])

# 2.1 for each cnv from cnvd, find the is around the boundary, and the same with tandom sites 
CNV_CNVD_InsulationScore_GM12878_HGPIN  <- numeric()
random_CNVD_insulationScore_GM12878_HGPIN <- numeric()
for ( i in names(table(CNV_CNVD[filter_id_3 & filter_id_4 & filter_id_6 , 'chr'])) ) {
  i <- as.numeric(i)
  tmp_sites <-  CNV_CNVD[filter_id_3 & filter_id_4 & filter_id_6 , 3:5]
  tmp_sites <- tmp_sites[tmp_sites[,1] == i, ]
  tmp_cnv_block_ins <- Sites_IS( Sites = tmp_sites, insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
  CNV_CNVD_InsulationScore_GM12878_HGPIN <- rbind(CNV_CNVD_InsulationScore_GM12878_HGPIN, tmp_cnv_block_ins)
  # tandom sites 
  tmp_random_number <- dim(tmp_sites)[1] # for each chromosome random select some sites with the same number of TADs
  tmp_random_site <- GM12878_is_10kb[[i]][ sample( x = 1:dim(GM12878_is_10kb[[i]])[1], size = tmp_random_number), ]
  tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
  random_CNVD_insulationScore_GM12878_HGPIN <- rbind(random_CNVD_insulationScore_GM12878_HGPIN,  tmp_random_ins)
}
median_ins_CNV_CNVD_HGPIN <- apply(CNV_CNVD_InsulationScore_GM12878_HGPIN, 2, function(x) median( as.numeric(x), na.rm = T) )
median_random_insulationScore_GM12878_CNVD_HGPINC <- apply(random_CNVD_insulationScore_GM12878_HGPIN, 2, function(x) median( as.numeric(x), na.rm = T) )

# 2.4 draw the insulation score around the cnv sites
png(filename = '../CNV_TAD_boundary/Human disease (HGPINC) cnv boundary and TAD boundary (GM12878).png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1, mar = c(3,3,3,2))
plot( median_ins_CNV_CNVD_HGPIN, type = 'l', 
      xlab = 'Genomic Regions', 
      ylab = 'Median Insulation Score',
      main = 'CNV Boundaries in Cancers \nand TAD boundaris in GM12878',
      ylim = c(min(median_ins_CNV_CNVD_HGPIN, GM12878_median_TAD_ins,median_random_insulationScore_GM12878_CNVD_HGPINC, na.rm = T), 
               max(median_ins_CNV_CNVD_HGPIN, GM12878_median_TAD_ins,median_random_insulationScore_GM12878_CNVD_HGPINC, na.rm = T) ),
      col = "red",
      lwd = 5, 
      #xlim = , 
      axes = FALSE)
lines( GM12878_median_TAD_ins, lwd = 5, col = "black" )
lines(median_random_insulationScore_GM12878_CNVD_HGPINC, lwd = 5, col = "blue")
axis(side = 1, at = c(1, 50, 101, 151, 201), labels = c('-1000kb', '-500kb','0kb','500kb', '1000kb' ))
axis(side = 2)
legend('bottomright', bty = 'n', cex = 0.6, 
       legend = c('CNV boundaries\n(HGPINC)', 'TAD Boudaries', 'Random Sites'),
       col = c( 'red', 'black', 'blue'), pch = 19)
dev.off()


# select a disease : Non-small cell lung cancer
filter_id_7 <- CNV_CNVD[, 'disease'] == 'Non-small cell lung cancer'
# CNV_CNVD[, 'disease'] == 'Basal cell lymphoma' | CNV_CNVD[, 'disease'] == 'Acute myeloid leukemia' | CNV_CNVD[, 'disease'] == 'Burkitts lymphoma' | 
sum(filter_id_7) # [1] 208
sum(filter_id_3 & filter_id_4 & filter_id_7)
head( CNV_CNVD[filter_id_3 & filter_id_4, ])

# 2.1 for each cnv from cnvd, find the is around the boundary, and the same with tandom sites 
CNV_CNVD_InsulationScore_GM12878_NSCLC  <- numeric()
random_CNVD_insulationScore_GM12878_NSCLC <- numeric()
for ( i in names(table(CNV_CNVD[filter_id_3 & filter_id_4 & filter_id_7 , 'chr'])) ) {
  i <- as.numeric(i)
  tmp_sites <-  CNV_CNVD[filter_id_3 & filter_id_4 & filter_id_7 , 3:5]
  tmp_sites <- tmp_sites[tmp_sites[,1] == i, ]
  tmp_cnv_block_ins <- Sites_IS( Sites = tmp_sites, insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
  CNV_CNVD_InsulationScore_GM12878_NSCLC <- rbind(CNV_CNVD_InsulationScore_GM12878_NSCLC, tmp_cnv_block_ins)
  # tandom sites 
  tmp_random_number <- dim(tmp_sites)[1] # for each chromosome random select some sites with the same number of TADs
  tmp_random_site <- GM12878_is_10kb[[i]][ sample( x = 1:dim(GM12878_is_10kb[[i]])[1], size = tmp_random_number), ]
  tmp_random_ins <- Sites_IS( Sites = tmp_random_site, insulation_score_list = GM12878_is_10kb[[i]], flank = 100)
  random_CNVD_insulationScore_GM12878_NSCLC <- rbind(random_CNVD_insulationScore_GM12878_NSCLC,  tmp_random_ins)
}
median_ins_CNV_CNVD_NSCLC <- apply(CNV_CNVD_InsulationScore_GM12878_NSCLC, 2, function(x) median( as.numeric(x), na.rm = T) )
median_random_insulationScore_GM12878_CNVD_NSCLC <- apply(random_CNVD_insulationScore_GM12878_NSCLC, 2, function(x) median( as.numeric(x), na.rm = T) )

# 2.4 draw the insulation score around the cnv sites
png(filename = '../CNV_TAD_boundary/Human disease (NSCLC) cnv boundary and TAD boundary (GM12878).png', width = 1024, height = 1024)
par(cex = 3, mgp = c(1.5, 0.5, 0), lwd = 3, tcl =-0.1, mar = c(3,3,3,2))
plot( median_ins_CNV_CNVD_NSCLC, type = 'l', 
      xlab = 'Genomic Regions', 
      ylab = 'Median Insulation Score',
      main = 'CNV Breakpoints in Cancers \nand TAD boundaris in GM12878',
      ylim = c(min(median_ins_CNV_CNVD_NSCLC, GM12878_median_TAD_ins,median_random_insulationScore_GM12878_CNVD_NSCLC, na.rm = T), 
               max(median_ins_CNV_CNVD_NSCLC, GM12878_median_TAD_ins,median_random_insulationScore_GM12878_CNVD_NSCLC, na.rm = T) ),
      col = "red",
      lwd = 5, 
      #xlim = , 
      axes = FALSE)
lines( GM12878_median_TAD_ins, lwd = 5, col = "black" )
lines(median_random_insulationScore_GM12878_CNVD_NSCLC, lwd = 5, col = "blue")
axis(side = 1, at = c(1, 50, 101, 151, 201), labels = c('-1000kb', '-500kb','0kb','500kb', '1000kb' ))
axis(side = 2)
legend('bottomright', bty = 'n', cex = 0.6, 
       legend = c('CNV Breakpoints\n(NSCLC)', 'TAD Boudaries', 'Random Sites'),
       col = c( 'red', 'black', 'blue'), pch = 19)
dev.off()

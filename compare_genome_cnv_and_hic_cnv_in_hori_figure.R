# compare genome cnv and hic cnv, 3 cell lines
# 20160201
# wupz
# library packages and source scripts
source('/lustre/user/liclab/wupz/dosageEffect/scripts/DrawCNVFigure.R')
# 1. read data
# 1.1 Data of RPMI-8226
RPMI8226_hic_cnv_raw <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/hic_cnv_raw/RPMI_8226_merged.bam_ratio.txt', header = T, stringsAsFactors = F)
RPMI8226_hic_cnv_filter <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/hic_cnv_filtered/RPMI8226_XY.bed_ratio.txt', header = T, stringsAsFactors = F)
RPMI8226_genome_cnv <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_200kb/8226-merged.sorted.bam.dedup.bam_ratio.txt',header = T,  stringsAsFactors = F)
dim(RPMI8226_hic_cnv_raw)
dim(RPMI8226_hic_cnv_filter)
dim(RPMI8226_genome_cnv)
# 1.2 Data of U266
U266_genome_cnv<- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_200kb/U266-merged.sorted.bam.dedup.bam_ratio.txt', header = T, stringsAsFactors = F)
U266_hic_cnv_raw <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/hic_cnv_raw/U266_merged.bam_ratio.txt', header = T, stringsAsFactors = F)
U266_hic_cnv_filter <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/hic_cnv_filtered/U266_XY.bed_ratio.txt', header = T, stringsAsFactors = F)

# 1.3 Data of GM12878
GM12878_genome_cnv<- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/GM12878_CNV/ERR091571.sam_ratio.txt', header = T, stringsAsFactors = F)
GM12878_hic_cnv_raw <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/GM12878_CNV/GM12878_hic_unfiltered_merge.bam_ratio.txt', header = T, stringsAsFactors = F)
GM12878_hic_cnv_filter <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/GM12878_CNV/GM12878_hic_filter.bed_ratio.txt', header = T, stringsAsFactors = F)

# 2 Draw the cnv data
# each cell line in a figure with 3 sub-figures
# figure dir: '/lustre/user/liclab/wupz/dosageEffect/hicCNV
# 2.1 RPMI-8226
png( filename = 'Comparing_the_HiC_CNV_and_genome_CNV_RPMI8226.png', width = 1536, height = 1280)
layout(mat = matrix(1:3, nrow = 3) )
par(cex = 3, lwd = 2, mar = c(1, 2, 2, 0.5), mgp = c(0.8,0,0), tcl = 0.1, cex.lab = 1, cex.axis = 0.8, font = 2, xaxs = 'i')
# RPMI8226, Genome CNV
DrawCNVFigure(CNV_data = RPMI8226_genome_cnv, ploidy = 3, mode = 'horizontal', ylab = '', main = '', x_axis_lab = F)
# RPMI8226, Raw Hi-C
par( mar = c(1, 2, 0.2, 0.5))
DrawCNVFigure(CNV_data = RPMI8226_hic_cnv_raw, ploidy = 3, mode = 'horizontal', ylab = '', x_axis_lab = F)
# RPMI8226, Filtered Hi-C
par(cex.axis = 1, font = 2)
DrawCNVFigure(CNV_data = RPMI8226_hic_cnv_filter, ploidy = 3, mode = 'horizontal', ylab = '', x_axis_lab = T, yaxt="n")
# add y axis
axis(2,cex.axis=1)
dev.off()


png( filename = 'Comparing_the_HiC_CNV_and_genome_CNV_U266.png', width = 1536, height = 1280)
layout(mat = matrix(1:3, nrow = 3) )
par(cex = 3, lwd = 2, mar = c(1, 2, 2, 0.5), mgp = c(0.8,0,0), tcl = 0.1, cex.lab = 1, cex.axis = 0.8, font = 2, xaxs = 'i')
# RPMI8226, Genome CNV
DrawCNVFigure(CNV_data = U266_genome_cnv, ploidy = 2, mode = 'horizontal', ylab = '',main = '', x_axis_lab = F)
# RPMI8226, Raw Hi-C
par( mar = c(1, 2, 0.2, 0.5))
DrawCNVFigure(CNV_data = U266_hic_cnv_raw, ploidy = 2, mode = 'horizontal', ylab = '', x_axis_lab = F)
# RPMI8226, Filtered Hi-C
par(cex.axis = 1, font = 2)
DrawCNVFigure(CNV_data = U266_hic_cnv_filter, ploidy = 2, mode = 'horizontal', ylab = '', x_axis_lab = T, yaxt="n")
# add y axis
axis(2,cex.axis=1)
dev.off()


png( filename = 'Comparing_the_HiC_CNV_and_genome_CNV_GM12878.png', width = 1536, height = 1280)
layout(mat = matrix(1:3, nrow = 3) )
par(cex = 3, lwd = 2, mar = c(1, 2, 2, 0.5), mgp = c(0.8,0,0), tcl = 0.1, cex.lab = 1, cex.axis = 0.8, font = 2, xaxs = 'i')
# RPMI8226, Genome CNV
DrawCNVFigure(CNV_data = GM12878_genome_cnv, ploidy = 2, mode = 'horizontal', ylab = 'WGS', main = 'Copy Number of GM12878', x_axis_lab = F)
# RPMI8226, Raw Hi-C
par( mar = c(1, 2, 0.2, 0.5))
DrawCNVFigure(CNV_data = GM12878_hic_cnv_raw, ploidy = 2, mode = 'horizontal', ylab = 'Raw Hi-C Data', x_axis_lab = F)
# RPMI8226, Filtered Hi-C
par(cex.axis = 0.6, font = 2)
DrawCNVFigure(CNV_data = GM12878_hic_cnv_filter, ploidy = 2, mode = 'horizontal', ylab = 'Filtered Hi-C Data', x_axis_lab = T, yaxt="n")
# add y axis
axis(2,cex.axis=0.8)
dev.off()


# 3. draw the smoothed scatter plot of genome cnv and filtered hic
# 3.1 RPMI-8226
png( filename = 'Scatter_plot_of_CNV_RPMI8226.png', width = 1024, height = 1024)
par(cex = 3, lwd = 2, mar = c(2, 2, 3.5, 0.5), mgp = c(0.8,0,0), tcl = 0.1, cex.lab = 1, cex.axis = 0.8, font = 2)
x_vec <- RPMI8226_genome_cnv[RPMI8226_genome_cnv$Ratio != -1 | RPMI8226_hic_cnv_filter$Ratio != -1,3]*3
y_vec <- RPMI8226_hic_cnv_filter[RPMI8226_genome_cnv$Ratio != -1 | RPMI8226_hic_cnv_filter$Ratio != -1, 3]*3
correlation_x_y <- cor(x_vec, y_vec, method = 'spearman')
cor.test(x_vec, y_vec, method = 'spearman') # p val: < 2.2 * 10^(-16)
tmp_pval <- expression( 2.2 * 10^(-16) )
smoothScatter( x = x_vec, 
               y = y_vec, 
               asp = 1,
               xlim = c(0,10),
               ylim = c(0,10),
               main= 'Copy Number of RPMI-8226 (200kb bins)\nWGS v.s. Filtered Hi-C Data',
               xlab = 'WGS Estimated Copy Numbers',
               ylab = 'Filtered Hi-C Estimated Copy Numbers'
               )
abline(a = 0, b = 1)
legend( 'topleft', legend = paste( 'Spearman Correlation:', round(correlation_x_y, 3), tmp_pval) , bty = 'n')
dev.off()
# 3.2 U266
png( filename = 'Scatter_plot_of_CNV_U266.png', width = 1024, height = 1024)
par(cex = 3, lwd = 2, mar = c(2, 2, 3.5, 0.5), mgp = c(0.8,0,0), tcl = 0.1, cex.lab = 1, cex.axis = 0.8, font = 2)
x_vec <- U266_genome_cnv[U266_genome_cnv$Ratio != -1 | U266_hic_cnv_filter$Ratio != -1,3]*2
y_vec <- U266_hic_cnv_filter[U266_genome_cnv$Ratio != -1 | U266_hic_cnv_filter$Ratio != -1, 3]*2
correlation_x_y <- cor(x_vec, y_vec, method = 'spearman')
cor.test(x_vec, y_vec, method = 'spearman') # p-value < 0.00000000000000022
smoothScatter( x = x_vec, 
               y = y_vec, 
               asp = 1,
               xlim = c(0,10),
               ylim = c(0,10),
               main= 'Copy Number of U266 (200kb bins)\nWGS v.s. Filtered Hi-C Data',
               xlab = 'Whole Genome Sequencing',
               ylab = 'Filtered Hi-C'
)
abline(a = 0, b = 1)
legend( 'topleft', legend = paste( 'Spearman Correlation:', round(correlation_x_y, 3)), bty = 'n')
dev.off()

# 3.3 GM12878
png( filename = 'Scatter_plot_of_CNV_GM12878.png', width = 1024, height = 1024)
par(cex = 3, lwd = 2, mar = c(2, 2, 3.5, 0.5), mgp = c(0.8,0,0), tcl = 0.1, cex.lab = 1, cex.axis = 0.8, font = 2)
x_vec <- GM12878_genome_cnv[GM12878_genome_cnv$Ratio != -1 | GM12878_hic_cnv_filter$Ratio != -1,3]*2
y_vec <- GM12878_hic_cnv_filter[GM12878_genome_cnv$Ratio != -1 | GM12878_hic_cnv_filter$Ratio != -1, 3]*2
correlation_x_y <- cor(x_vec, y_vec, method = 'spearman')
smoothScatter( x = x_vec, 
               y = y_vec, 
               asp = 1,
               xlim = c(0,10),
               ylim = c(0,10),
               main= 'Copy Number of GM12878 (200kb bins)\nWGS v.s. Filtered Hi-C Data',
               xlab = 'Whole Genome Sequencing',
               ylab = 'Filtered Hi-C'
)
legend( 'topleft', legend = paste( 'Spearman Correlation:', round(correlation_x_y, 3)), bty = 'n')
abline(a = 0, b = 1)
dev.off()

# END
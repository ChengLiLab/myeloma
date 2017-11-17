# compare genome cnv and hic cnv, plot as a 5*5 matrix
# 20151229
# wupz
# library packages and source scripts
source('/lustre/user/liclab/wupz/dosageEffect/scripts/DrawCNVFigure.R')
# 1. read data
RMPI8226_hic_cnv_raw <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/hic_cnv_raw/RMPI_8226_merged.bam_ratio.txt', header = T)
RMPI8226_hic_cnv_filter <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/hic_cnv_filtered/RMPI8226_XY.bed_ratio.txt', header = T)
RMPI8226_genome_cnv <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_200kb/8226-merged.sorted.bam.dedup.bam_ratio.txt', 
                                  header = T,  stringsAsFactors = F)
dim(RMPI8226_hic_cnv_raw)
dim(RMPI8226_hic_cnv_filter)
dim(RMPI8226_genome_cnv)
U266_genome_cnv<- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_200kb/U266-merged.sorted.bam.dedup.bam_ratio.txt', header = T)
U266_hic_cnv_raw <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/hic_cnv_raw/U266_merged.bam_ratio.txt', header = T)
U266_hic_cnv_filter <- read.table('/lustre/user/liclab/wupz/DNA_seq/CNV_final/hic_cnv_filtered/U266_XY.bed_ratio.txt', header = T)

# 2. plot the genome cnv and hic cnv (raw ) in a single figure
DrawCNVFigure(CNV_data = RMPI8226_hic_cnv_raw, ploidy = 3, title = 'RMPI8226')
# compare two cnv dataset

Draw2CNVFigure( CNV_data_1 = RMPI8226_genome_cnv, 
                CNV_data_2 = RMPI8226_hic_cnv_filter, 
                ploidy = c(3,3), 
                legends = c('Genome-Seq', 'Filter Hi-C'), 
                title= 'RMPI-8226', 
                filename = 'Compare_genome_cnv_and_filter_hic_cnv.png')
Draw2CNVFigure( CNV_data_1 = RMPI8226_genome_cnv, 
                CNV_data_2 = RMPI8226_hic_cnv_raw, 
                ploidy = c(3,3), 
                legends = c('Genome-Seq', 'Raw Hi-C'), 
                title= 'RMPI-8226', 
                filename = 'Compare_genome_cnv_and_raw_hic_cnv.png')

Draw2CNVFigure( CNV_data_1 = U266_genome_cnv, 
                CNV_data_2 = U266_hic_cnv_filter, 
                ploidy = c(3,3), 
                legends = c('Genome-Seq', 'Filter Hi-C'), 
                title= 'U266', 
                filename = 'Compare_genome_cnv_and_filter_hic_cnv_U266.png')
Draw2CNVFigure( CNV_data_1 = U266_genome_cnv, 
                CNV_data_2 = U266_hic_cnv_raw, 
                ploidy = c(2,2), 
                legends = c('Genome-Seq', 'Raw Hi-C'), 
                title= 'U266', 
                filename = 'Compare_genome_cnv_and_raw_hic_cnv_U266.png')
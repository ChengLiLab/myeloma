# test PlotTwoChromosomeHic
# 20161123
# wupz
PlotTwoChromosomeHic( two_hic_matrix = RMPI_8226_HindIII_chr16_chr22, 
                      two_cnv = RPMI8226_CNV_chr16_chr22, 
                      ploidy = 3, 
                      chromosome_name = c("chr16", "chr22"), 
                      chromosome_length = c(hg19_len[16, 3], hg19_len[22, 3]),
                      cex = 1 )
dev.off()
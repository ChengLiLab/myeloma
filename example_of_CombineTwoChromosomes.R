# test CombineTwoChromosomes
# 20161123
# wupz
#chr5_chr10 <- CombineTwoChromosomes( raw_matrix = raw_HiC_matrix_2M[[1]], 5,10)

PlotHicHeatmap <- function( matrix, title = '', legend_lab = '') {
  Heatmap(matrix, 
          name = legend_lab, 
          column_title = title,
          column_title_gp = gpar(fontsize = 100), 
          col = colorRamp2( quantile(matrix, c(0.05, 0.95), na.rm = T), c( "white", "red") ),
          cluster_rows = F,
          cluster_columns = F,
          show_row_names = F,
          show_column_names = F,
          heatmap_legend_param = list(title_gp = gpar(fontsize = 75), labels_gp = gpar(fontsize = 75), color_bar = "continuous")
  )
  
}
# test
PlotHicHeatmap(chr5_chr10[1:10,1:10] )
# plot trans Hi-C matrix at possible site of translocation
# 20161010
# wupz

PlotTransHicMatrix <- function (trans_hic_matrix, chr_label1, chr_label2,
                          position1, position2,
                          flank = 1000000, resolution = 40000,
                          figure_name = NULL, ...) {
  if ( !is.null(figure_name) ) {
    png( filename = figure_name, width = 1024*7/6, height = 1024 )
    # set layout
    layout(mat = matrix(1:2, nrow = 1), widths = c(6, 1) )
  }
  else {
    layout(mat = matrix(1:2, nrow = 1), widths = c(6, 1) )
  }
  
  # 1.1 plot trans_hic_matrix
  par(mar = c(3,3,3,3))
  position1 <- as.integer(position1)
  position2 <- as.integer(position2)
  chr1_start <- floor(position1/resolution) - flank/resolution
  chr1_end <- floor(position1/resolution) + flank/resolution
  chr2_start <- floor(position2/resolution) - flank/resolution
  chr2_end <- floor(position2/resolution) + flank/resolution
  trans_hic_matrix_part <- trans_hic_matrix[chr1_start:chr1_end, chr2_start:chr2_end]
  print(c(chr1_start, chr1_end, chr2_start, chr2_end) )
  # 1.2 plot trans_hic_matrix
  print(par()$mar)
  rect_heatmap(trans_hic_matrix_part, xlab = "", ylab = "", ... ) 

  # add lines
  end_line_pos <- flank/resolution * 2 + 1 
  center_pos <- flank/resolution + 1
  segments( x0 = c(0, end_line_pos), y0 = 0, x1 = c(0, end_line_pos), y1 = -end_line_pos, lty = 2, ...)
  segments( x0 = 0, y0 =  c(0, -end_line_pos), x1 = end_line_pos, y1 =  c(0, -end_line_pos), lty = 2, ...)
  # add labels of the vertical chromosomes
  mtext( chr_label1, side = 2, at = -center_pos, line = -6, ... )
  mtext( paste( round(c( position1 - flank, position1, position1 + flank)/1000000, 2), "M", sep = ""), 
         side = 2, at = c(0, -center_pos, -end_line_pos), 
         line = -1.5, adj = c(1, Inf, 0), ...  ) # use Inf to indicate center alignment
  # add labels of the horizontal chromosomes
  mtext( chr_label2, side = 1, at = center_pos, line = - 5, ... )
  mtext( paste( round(c( position2 - flank, position2, position2 + flank)/1000000, 2), "M", sep = ""), 
         side = 1, at = c(0, center_pos, end_line_pos), 
         line = 0, adj = c(0, Inf, 1), ... )
  # 2. plot heatmap legend
  par( mar = c(5, 0, 5, 1) )
  add_legend( trans_hic_matrix_part, ...)
  # 3. save the figure
  if ( !is.null(figure_name) ) {
    dev.off()
  }
}


# test PlotTransHicMatrix 
PlotTransHicMatrix(trans_hic_matrix = RMPI_8226_HindIII_chr16_chr22_40kb, 
                   chr_label1 = "chr16", chr_label2 = "chr22", 
                   position1 = 79966531, position2 = 23264891, 
                   flank = 1000000, resolution = 40000, cex = 3, lwd = 10, 
                   figure_name = "testPlotTransHicMatrix.png")





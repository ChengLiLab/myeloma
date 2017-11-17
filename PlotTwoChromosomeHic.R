# Function script: PlotTwoChromosomeHic
# plot Hi-C matrix of two chromosomes with CNV, 
# to check the translocation between the two chromosome
# 20161009
# wupz

# 1. source scripts
# setwd("/lustre/user/liclab/wupz/dosageEffect/hic_translocation")
source( "../scripts/CombineTwoChromosomes.R")
source( "../scripts/rect_heatmap.R")

# 2. set plot parameters
PlotTwoChromosomeHic <- function(two_hic_matrix, two_cnv, ploidy = 2, 
                                 chromosome_name = NULL, 
                                 chromosome_length = NULL,
                                 resolution = 200000,
                                 figure_name = NULL, ... ) {
  if ( !is.null(figure_name) ) {
    png( filename = figure_name, width = 2048 * 7/6, height = 2048 * 5/4)
    layout(mat = matrix( c(3,4,1,2), nrow = 2), widths = c(1,6), heights = c(4, 1))
  }
  # set layout, 1: two hic matrix; 2: cnv; 3: heatmap legend; 4: null
  layout(mat = matrix( c(3,4,1,2), nrow = 2), widths = c(1,6), heights = c(4, 1))
  par(mar = c(0,6,0,2) )
  # 1. plot Hi-C matrix of two chromosome
  rect_heatmap(raw_matrix = two_hic_matrix, 
               yaxt="n",ylab = "", xaxt="n", xlab = "", ...)
  # plot dashes lines to separate two chromosome
  print("Add separate line")
  dim_1 <- ceiling( chromosome_length[1]/resolution )
  dim_2 <- dim(two_hic_matrix)[1]
  # add verticality dash lines
  segments( x0 = c(0,  dim_2), y0 = c(0, -dim_1), x1 = c(0, dim_2), y1 = c(-dim_1, -dim_2), lty = 2, ... )
  # add horizontal dash line
  segments( x0 = c(0, dim_1), y0 =  c(0, -dim_2), x1 = c(dim_1, dim_2), y1 =  c(0, -dim_2), lty = 2, ... )
  # add blue dashed line to the trans-hic
  # # add verticality dash lines, blue
  segments( x0 = c(0, dim_1, dim_2), y0 = c(-dim_1, 0, 0), x1 = c(0, dim_1, dim_2), y1 = c(-dim_2, -dim_2, -dim_1), lty = 2,col = "blue", ...)
  # add horizontal dash line, blue
  segments( x0 = c(dim_1, 0, 0), y0 = c(0, -dim_1, -dim_2), x1 = c(dim_2, dim_2, dim_1), y1 = c(0, -dim_1, -dim_2), lty = 2, col = "blue", ...)
  # add y-axis labels corresponds to the chromosome length
  mtext( text = c(chromosome_name[1], chromosome_name[2]), 
         side = 2, at = c( -dim_1/2, -dim_1/2 - dim_2/2 ),
         line = -0.5,
         ... )
  mtext(text = c( "0M",  "0M"),
        side = 2, 
        at = c(0, -dim_1 - 5), 
        adj = 1, line = -0.5,
        ...
        )
  mtext(text = paste( round( c(chromosome_length[1], chromosome_length[2])/1000000 ), "M", sep = ""),
        side = 2, 
        at = c( - dim_1 + 5, -dim_2), 
        adj = 0, line = -0.5,
        ...
  )

  # add compartment A/B
  #par(mar = c(6,5,0,0), mgp = c(1,0.6,0), tcl = 1 )
  #plot(x = 0:(dim(RMPI_8226_HindIII_chr16_chr22_pca)[1]-1) * 2.5, 
  #     y = as.numeric( RMPI_8226_HindIII_chr16_chr22_pca[[1]] ), 
  #     type = "h", bty = "n",yaxt="n",ylab = "", xaxt="n", xlab = "",lwd = 10,
  #     col = ifelse(  as.numeric( RMPI_8226_HindIII_chr16_chr22_pca[[1]] ) >=0, "red", "blue") )
  #axis( side = 2, at = c( -0.05, 0.05), labels = c( "B", "A"), cex.axis = 6)
  
  # 2. plot cnv of two chromosomes
  par(mar = c(8,6,1,2), mgp = c(6, 3, 0))
  plot(x = 0:(dim_2-1) + 1/2, y = two_cnv[, 3] * ploidy, 
       ylim = c(0,6), pch = 19, bty = "n", yaxt="n", ylab = "", xaxt="n", xlab = "", 
       col = ifelse( two_cnv[, 5] == 2, "green", ifelse( two_cnv[, 5] > 2, "red", "blue")) )
  segments( x0 = dim_1, y0 =  0, x1 = dim_1, y1 =  6, lty = 2, ... )
  # add y-axis labels
  axis( side = 2, labels = F, line = -2, lwd = 10, lwd.ticks = 10)
  mtext( text = c(0, 2, 4, 6), at = c(0, 2, 4, 6), side = 2, line = -0.5, ...)
  # add x-axis labels
  axis( side = 1, at = c(0, dim_1/2, dim_1, dim_1/2 +dim_2/2, dim_2), labels = F, lwd = 10, lwd.ticks = 10 )
  mtext( text = c(chromosome_name[1], chromosome_name[2]), 
         side = 1, at = c(dim_1/2, dim_1/2 +dim_2/2), line = 5, 
         ... )
  mtext(text = c( "0 M",  "0M"),
        side = 1, 
        at = c(0, dim_1+5), 
        adj = 0, line = 5,
        ...
  )
  print( paste( round( c(chromosome_length[1], chromosome_length[2])/1000000 ), "M", sep = "") )
  mtext(text = paste( round( c(chromosome_length[1], chromosome_length[2])/1000000 ), "M", sep = ""),
        side = 1, 
        at = c(dim_1-5, dim_2), 
        adj = 1, line = 5,
        ...
  )
  # 3. add color legend of heatmap
  print("Add heatmap legend")
  par(mar = c(1,1,1,1))
  add_legend(two_hic_matrix, ...)
  if ( !is.null(figure_name) ) {
    dev.off()
  }
  print("Done")
}


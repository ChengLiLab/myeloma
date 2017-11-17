# A function to draw cnv data which generated from freec
# modified from freec
# 20151223
# wupz
DrawCNVFigure <- function (CNV_data , ploidy = 2, mode = 'squared', lines_h = T, lines_v = T,  filename = NULL,legend = '',x_axis_lab = T, ...) { # ylab = '', xlab = '',legend = '',main = '', 
  # input:
  #   CNV_data: the results of freec
  # ouput:
  #   a plot of CNV ratio and estimated copy number
  # para:
  #   mode: 
  #     'squared' or 'horizontal'. the former is a figure with 5 * 5 sub-figures each of a chromosomes, 
  #     the 'horizontal' mode is a plot with each chromosomes side by side
  #   filename: set the file names for saving the figure
  #   ploidy: set the ploidy of the normal genome
  #   main_lab: set the figure title
  ploidy <- ploidy
  #CNV_data$CopyNumber[CNV_data$Ratio == -1] <- 0
  if ( mode == 'squared' ) { # draw cnv data in a figure with 5 * 5 sub-figures
    if ( !is.null(filename) ) {
      png(filename = filename, width = 1280, height = 1280, units = "px", pointsize = 20, bg = "white", res = NA)
    }
    par(mfrow = c(5,5), mar = c( 1,1,1,1), 
        cex = 1.5, cex.main = 0.6, cex.lab = 0.6, cex.axis = 0.7, # set point size and labels fontsize
        mgp = c(0.5, 0.1, 0),
        tcl =-0.1
    )
    plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
    legend('center', bty = 'n', legend = paste('CN', c('<','=','>'),ploidy ), main = main,
           cex = 0.55, col = c('blue', 'green', 'red'), pch = 19)
    for (i in c(1:22,'X','Y')) {
      chr_i <- CNV_data$Chromosome == i
      chr_i_data <- as.data.frame(CNV_data[chr_i, ] )
      if ( length(chr_i ) > 0 ) {
        cnv_color <- ifelse( CNV_data[,5] ==2, 'green', ifelse(CNV_data[,5] > 2, 'red', 'blue') )
        plot( x = chr_i_data$Start, 
              y = chr_i_data$Ratio * ploidy,
              ylim = c(0, 6),
              main = paste ("Chromosome",i), 
              pch = '.',
              xaxt = "n", # rename xlab
              col = cnv_color,
              xlab = '',
              ylab = '')
        axis( side = 1, 
              at = round(quantile( range(chr_i_data$Start ), c(0, 0.5, 1) )),
              labels = paste( round( quantile( range(chr_i_data$Start ), c(0, 0.5, 1) ) / 1000000), 'M', sep = ' ' ) )
      }
    }
  }
  else if ( mode == 'horizontal') {
    if ( !is.null(filename) ) {
      png(filename = filename, width = 1280, height = 320, units = "px", pointsize = 20, bg = "white", res = NA)
      par(xaxs = 'i')
    }
    # draw the first chromosome cnv data: chr 1
    # set the color of cnv: (cn < 2, blue), (cn == 2, green), (cn > 2, red)
    col_ind <- ifelse( CNV_data[,5] ==2, 'green', ifelse(CNV_data[,5] > 2, 'red', 'blue') )
    col_cn_ind <- ifelse(CNV_data[,5] == 0, 'purple', 'black')
    # set the start and end coordinate of the chr1
    chr_start <- 1
    chr_end <- sum( CNV_data[,1] == 1)
    col_axis <- rep(c('black','red'),12)
    names(col_axis) <- c(1:22,'X', 'Y')
    plot(CNV_data[CNV_data[,1] == 1, 3] * ploidy, pch = '.', col = col_ind[chr_start:chr_end],
         xlim = c(0, dim(CNV_data)[1]),  ylim = c(0, 6),
         xaxt = "n", # rename axis
         ...
         )
    # set horizontal lines at each copy number
    if ( lines_h == T ) {
      abline( h = 0:6, lty = 3, col = 'grey30')
    }
    # set the start AND END boundary of chr1
    abline( v = chr_start, lty = 3, col = 'grey20')
    abline( v = chr_end, lty = 3, col = 'grey20')
    # Plot the estimated copy number
    points( x = chr_start:chr_end , y = CNV_data[CNV_data[,1] == 1, 5], pch = '.', col = col_cn_ind[CNV_data[,1] == 1] )
    # if x_axis_lab is true then plot the lab of x axis
    if( x_axis_lab ) {
      axis( side = 1, at = (chr_start + chr_end)/2, labels = 'Chr 1', line=0, col.axis = col_axis[1])
    }
    
    legend('topleft', legend = legend, bty = 'n')

    for( i in c(2:22,'X', 'Y')) {
      if(sum( CNV_data[,1] == i) != 0) {
        chr_start <- chr_end + 1
        chr_end <- chr_end + sum( CNV_data[,1] == i)
        abline( v = chr_end, lty = 3, col = 'grey20')
        points( x = chr_start:chr_end , y = CNV_data[CNV_data[,1] == i, 3] * ploidy, pch = '.',col = col_ind[CNV_data[,1] == i])
        points( x = chr_start:chr_end , y = CNV_data[CNV_data[,1] == i, 5], pch = '.', col = col_cn_ind[CNV_data[,1] == i] )
        # add x axis lab
        if( x_axis_lab & is.element(i, c(seq(3, 21, 2), "X") ) ) {
          axis( side = 1, at = (chr_start + chr_end)/2, labels = i, line=0, col.axis = col_axis[i] )
        }
      }
    }
  }
  else { 
    print('mode must be one of "squared" or "horizontal"')
    }
  if ( !is.null(filename) ) {
    dev.off()
  }
}
#test
#DrawCNVFigure(CNV_data = RMPI8226_hic_cnv_raw, ploidy = 3,ylab = '', mode = 'horizontal', x_axis_lab = T, lines_h = F)

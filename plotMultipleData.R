# Plot Hi-C, ChIP-seq, RNA-seq data
# 20160822
# wupz
plotMultipleData <- function( Hic_matrix, 
                              Chipseq_list, 
                              RNAseq_list = NULL, 
                              ylabs, 
                              main = NULL,
                              filename = NULL,
                              genomic_range = NULL) {
  source("../scripts/tri_heatmap.r")
  
  Chipseq_length <- length(Chipseq_list)
  RNAseq_length <- length(RNAseq_list)
  # dev.new()
  if ( !is.null(filename) ) {
    png(filename = filename, width = 1538, height = 2048)
    par(cex = 3)
  }
  layout(mat = matrix( c( 1:( 1 + Chipseq_length + RNAseq_length) ), 
                       nrow = 1 + Chipseq_length + RNAseq_length ), 
         heights = c(5, rep.int(1, Chipseq_length + RNAseq_length))  )
  # plot Hi-C matrix
  par(mar = c(1, 3, 2, 1), cex = 2)
  tri_heatmap(Hic_matrix, main = main)
  # set ChIP-seq colors
  palette(rainbow(10))
  # plot ChIP-seq peaks, format: bed
  
  for ( i in 1:Chipseq_length ) {
    par(mar = c(1, 3, 0, 1), mgp = c(1,0.5,0), cex = 2 ) 
    if( is.null(genomic_range) ) {
      xlim_1 <- c(0, max(Chipseq_list[[i]][, 3]) )
    }
    else if ( !is.null(genomic_range) ) {
      xlim_1 <- genomic_range
    }
    plot(x = xlim_1, 
         y = c(0, max( log2(Chipseq_list[[i]][, 5] ), na.rm = T)*1.05 ), 
         type = "n", xlab = "", xaxt = "n", # suppress x axes and lab
         ylab = ylabs[i], yaxt = "n", bty = "n",
         xlim = xlim_1
         )
    rect(xleft = Chipseq_list[[i]][, 2], 
         ybottom = 0,
         xright = Chipseq_list[[i]][, 3],
         ytop = log2( Chipseq_list[[i]][, 5] ),
         col = palette()[i], border = palette()[i]
         )
  }
  # plot RNA-seq peaks, format: bed
  if( !is.null(RNAseq_list) ) {
    for ( i in 1:RNAseq_length ) {
      par(mar = c(1, 3, 0, 1), mgp = c(1,0.5,0),cex = 2 ) 
      if( is.null(genomic_range) ) {
        xlim_1 <- c(0, max(RNAseq_list[[i]][, 3]) )
      }
      else if ( !is.null(genomic_range) ) {
        xlim_1 <- genomic_range
      }
      plot(x = xlim_1, 
           y = c(0, max( log2(RNAseq_list[[i]][, 5]), na.rm = T )*1.05 ), 
           type = "n", xlab = "", xaxt = "n", # suppress x axes and lab
           ylab = ylabs[i + Chipseq_length ], yaxt = "n", bty = "n") 
      rect( xleft = RNAseq_list[[i]][, 2], 
            ybottom = 0,
            xright = RNAseq_list[[i]][, 3],
            ytop = log2( RNAseq_list[[i]][, 5] ),
            col = palette()[i], border = palette()[i]
            )
    }
  }
  if ( !is.null(filename) ) {
    dev.off()
  }
}
# Plot Hi-C, ChIP-seq, RNA-seq data
# 20160822
# wupz
plotMultipleData2 <- function( Hic_matrix_list, 
                              TAD_insulation_score_list = NULL,
                              freec_CNV = NULL, ploidy = 2,
                              Chipseq_list = NULL, 
                              RNAseq_list = NULL, 
                              ylabs, 
                              main = NULL,
                              filename = NULL,
                              genomic_range = NULL,
                              resolution = 200000,
                              gene_legend = c(name = NULL, posi = NULL), ... ) {
  library(limma)
  source("../scripts/tri_heatmap.r")
  if( is.null(TAD_insulation_score_list) ) {
    plot_TAD <- 0
  } else {
    plot_TAD <- 1
  }
  if( is.null(freec_CNV) ) {
    plot_cnv <- 0
  } else {
    plot_cnv <- 1
  }
  Chipseq_length <- length(Chipseq_list)
  RNAseq_length <- length(RNAseq_list)
  # dev.new()
  if ( !is.null(filename) ) {
    png(filename = filename, width = 1538, height = 2048)
    par(... )
  }
  # check genome range
  if( is.null(genomic_range) ) {
    print("genomic_range is null! ")
  }
  else if ( !is.null(genomic_range) ) {
    tmp_chr <- strsplit2(genomic_range, split = ":|-")[1]
    tmp_chr_start <- as.numeric( strsplit2(genomic_range, split = ":|-")[2] ) + 1
    tmp_chr_end <- as.numeric( strsplit2( genomic_range, split = ":|-")[3] )
    # calculate the chromosome id and ranges
    if ( tmp_chr == "chrX") {
      tmp_chr_id = 23
    }
    else if ( tmp_chr == "chrY") {
      tmp_chr_id = 24
    }
    else {
      tmp_chr_id <- as.integer( strsplit2( x = tmp_chr, split = "chr")[2] )
    }
    tmp_chr_start_id <- ceiling( tmp_chr_start/resolution ) 
    tmp_chr_end_id <- ceiling( tmp_chr_end/resolution ) 
    # fix the xlim range
    tmp_xlim <- c(tmp_chr_start_id - 1, tmp_chr_end_id )
    # print(c(tmp_chr, tmp_chr_start,tmp_chr_end, tmp_xlim))
  }
  layout(mat = matrix( c( 1:( 1 + plot_TAD + plot_cnv + Chipseq_length + RNAseq_length) ), 
                       nrow = 1 + plot_TAD + plot_cnv + Chipseq_length + RNAseq_length ), 
         heights = c(3, rep.int(1, plot_TAD + plot_cnv + Chipseq_length + RNAseq_length))  )
  # plot Hi-C matrix
  par(mar = c(1, 2, 2, 1), mgp = c(1.5,0.5,0), ... )
  tmp_hic_matrix <- Hic_matrix_list[[tmp_chr_id]][[tmp_chr_id]][tmp_chr_start_id:tmp_chr_end_id, tmp_chr_start_id:tmp_chr_end_id]
  tri_heatmap(raw_matrix = tmp_hic_matrix, main = main, cex.main = 2)
  # set ChIP-seq colors
  palette(rainbow(7))
  # plot TAD and insulation score 
  if ( plot_TAD ) {
    par(mar = c(1, 2, 0, 1), mgp = c(1.5,0.5,0), ...)
    # find insulation score data
    tmp_insulation_score_data_x <- TAD_insulation_score_list[[1]][tmp_chr_start_id:tmp_chr_end_id, 4]
    tmp_insulation_score_data_y <- TAD_insulation_score_list[[1]][tmp_chr_start_id:tmp_chr_end_id, 8]
    tmp_TAD_insulation_score_data_id <- (TAD_insulation_score_list[[2]][, 2] >= tmp_chr_start) & 
      (TAD_insulation_score_list[[2]][, 3] <= tmp_chr_end)
    tmp_TAD_insulation_score_data_x <- TAD_insulation_score_list[[2]][tmp_TAD_insulation_score_data_id, 2] + resolution/2
    tmp_TAD_insulation_score_data_y <- TAD_insulation_score_list[[2]][tmp_TAD_insulation_score_data_id, 8] 
    # plot
    plot(x = tmp_insulation_score_data_x, 
         y = tmp_insulation_score_data_y,
         xlim = tmp_xlim*resolution,
         ylim = c(-2,2), 
         cex = 3, col = "black", type = "l", lwd = 8,
         xlab = "", xaxt = "n",yaxt = "n", bty = "n", ylab = ""
    )
    segments(x0 = tmp_TAD_insulation_score_data_x, x1 = tmp_TAD_insulation_score_data_x, 
             y0 = -2, y1 = 2, 
             lwd = 8, col= "red")
    axis(side = 2, at = c(-2,0,2), labels = c(-2,0,2) )
  }

  # plot cnv data
  if ( plot_cnv ) {
    par(mar = c(1, 2, 0, 1), mgp = c(1.5,0.5,0), ... )
    tmp_cnv_x <- freec_CNV[freec_CNV[, 1] == tmp_chr_id, 2][tmp_chr_start_id:tmp_chr_end_id] + resolution/2
    tmp_cnv_y <- freec_CNV[freec_CNV[, 1] == tmp_chr_id, 3][tmp_chr_start_id:tmp_chr_end_id] * ploidy
    tmp_cnv_y_col <- freec_CNV[freec_CNV[, 1] == tmp_chr_id, 5][tmp_chr_start_id:tmp_chr_end_id] 
    plot(x = tmp_cnv_x, 
         y = tmp_cnv_y,
         xlim = tmp_xlim*resolution,
         ylim = c(0,6), pch = 16,
         col = ifelse(tmp_cnv_y_col ==2, "green", ifelse(tmp_cnv_y_col>2, "red", "blue")), 
         xlab = "", xaxt = "n",yaxt = "n", bty = "n", ylab = ""
    )
    axis(side = 2, at = c(0,2,4,6), labels = c(0,2,4,6), cex.axis= 0.8)
  }
  
  # plot ChIP-seq peaks, format: bed
  if ( Chipseq_length ) {
    for ( i in 1:Chipseq_length ) {
      par(mar = c(1, 2, 0, 1), mgp = c(1.5,0.5,0), ... )
      # found chiq-seq data
      tmp_chipseq_chr <- Chipseq_list[[i]][, 1] == tmp_chr
      
      # print( head(Chipseq_list[[i]][tmp_chipseq_chr, ]))
      tmp_chipseq_id <- (Chipseq_list[[i]][tmp_chipseq_chr, 2] >= tmp_chr_start) & 
        (Chipseq_list[[i]][tmp_chipseq_chr, 3] <= tmp_chr_end)
      
      tmp_chipseq <- Chipseq_list[[i]][tmp_chipseq_chr, ][tmp_chipseq_id, ]
      if ( !is.infinite(max( log2(tmp_chipseq[, 5] ), na.rm = T)) ) {
        tmp_ylim <- c(0, max( log2(tmp_chipseq[, 5] ), na.rm = T)*1.05 )
      } else {
        tmp_ylim <- c(0, 10)
      }
      
      plot(x = tmp_xlim, 
           y = tmp_ylim, 
           type = "n", xlab = "", xaxt = "n", # suppress x axes and lab
           ylab = ylabs[i], yaxt = "n", bty = "n"
      )
      if( sum(tmp_chipseq_id) > 0 ) {
        rect(xleft = tmp_chipseq[, 2]/resolution, 
             ybottom = 0,
             xright = tmp_chipseq[, 3]/resolution,
             ytop = log2( tmp_chipseq[, 5] ),
             col = palette()[i], border = palette()[i]
        )
      }
      
      # add y axis and ylab
      print(floor(tmp_ylim * 0.8))
      axis(side = 2, at = round(tmp_ylim, digits = 1 ), labels = round(tmp_ylim, digits = 1 ) )
      # test borders
      # abline(v = tmp_xlim, col = "black", lwd = 3)
      # abline(v =  min(tmp_chipseq[, 2]/resolution), col = "blue", lwd = 3)
      # print( min(tmp_chipseq[, 2]/resolution))
      # abline(v =  max(tmp_chipseq[, 3]/resolution), col = "blue", lwd = 3)
      # print(max(tmp_chipseq[, 3]))
    }
  }
  
  # plot RNA-seq peaks, format: bed
  if( !is.null(RNAseq_list) ) {
    for ( i in 1:RNAseq_length ) {
      par(mar = c(2, 2, 0, 1), mgp = c(1.5,0.5,0), ... )
      # found RNA-seq data
      tmp_RNAseq_chr <- RNAseq_list[[i]][, 1] == tmp_chr
      tmp_RNAseq_id <- (RNAseq_list[[i]][tmp_RNAseq_chr, 2] >= tmp_chr_start) & 
        (RNAseq_list[[i]][tmp_RNAseq_chr, 3] <= tmp_chr_end)
      tmp_RNAseq <- RNAseq_list[[i]][tmp_RNAseq_chr, ][tmp_RNAseq_id, ]
      tmp_ylim <- c(0, 5)
      plot(x = tmp_xlim, 
           y = tmp_ylim, 
           type = "n", xlab = "", xaxt = "n", # suppress x axes and lab
           ylab = ylabs[i + Chipseq_length ], yaxt = "n", bty = "n") 
      rect( xleft = tmp_RNAseq[, 2]/resolution, 
            ybottom = 0,
            xright = tmp_RNAseq[, 3]/resolution,
            ytop = log2( tmp_RNAseq[, 5] ),
            col = palette()[i], border = palette()[i]
            )
    }
    print( min(tmp_RNAseq[, 2]/resolution))
    print(max(tmp_RNAseq[, 3]/resolution))
    axis(side = 2, at = tmp_ylim, labels = tmp_ylim )
  }
  if ( !is.null(gene_legend) ) {
    axis(side = 1, at = as.numeric( gene_legend[2] )/resolution, labels = gene_legend[1])
  }
  
  if ( !is.null(filename) ) {
    dev.off()
  }
}
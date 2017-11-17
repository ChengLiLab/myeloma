PlotBiasOfCNV <- function( cnv_block, cnv_block_mean,log2scale = T, legend_site = "bottom", file_name = NULL, ...) {
  # delete zeros and NAs
  tmp_non_zeros <- cnv_block_mean !=0 & !is.na(cnv_block_mean )
  # only show the copy number of 1 ~ 6
  tmp_cnv_1_6_ids <- cnv_block[, 5] >=1 &cnv_block[, 5] <= 6
  if( !is.null(file_name) ) {
    png(filename = file_name, width = 1024, height = 1024)
    par( cex = 3, mgp = c(1.5, 0.5, 0), pch = '.', lwd = 3, mar = c(3, 4, 4, 3))
  }
  tmp_y <- cnv_block_mean[tmp_non_zeros & tmp_cnv_1_6_ids] 
  tmp_x <- cnv_block[tmp_non_zeros & tmp_cnv_1_6_ids, 5]
  if( log2scale) {
    boxplot(log2(tmp_y) ~ tmp_x, col = c('blue', 'green', 'red', 'red', 'red', 'red'), 
            xlab = 'Copy Number of CNV Blocks',
            #ylab = 'Mean Hi-C Interaction Score\n log2 Scale',
            ... )
    # add linear regression line
    
    tmp_lm <- summary( lm( log2(tmp_y) ~ tmp_x ) )
  }
  else if ( !log2scale) {
    boxplot( tmp_y ~ tmp_x, col = c('blue', 'green', 'red', 'red', 'red', 'red'), 
             xlab = 'Copy Number of CNV Blocks',
             #ylab = 'Mean Hi-C Interaction Score\n log2 Scale',
             ... )
    # add linear regression line
    tmp_lm <- summary( lm( tmp_y ~ tmp_x ) )
  }
  tmp_slope <- coef(tmp_lm)[2, 1]
  tmp_intercept <- coef(tmp_lm)[1, 1]
  abline( a = tmp_intercept, b = tmp_slope, col = 'black' )
  legend(legend_site, bty = 'n', cex = 0.8,
         legend = paste('Slope of fitted linear regression line: ', round(coef(tmp_lm)[2, 1], 4), '\nP value: ', round(coef(tmp_lm)[2, 4], 10) ) 
  )
  if( !is.null(file_name) ) {
    dev.off()
  }
}
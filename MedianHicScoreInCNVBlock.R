# Calculate median hic score of a set of cnv block
# wupz
# 20160127
MedianHicScoreInCNVBlock <- function ( cnv_block, hic_matrix_list, methods = 'median') {
  hic_matrix_of_cnv_block <- list()
  for ( i in 1:dim(cnv_block)[1]) {
    tmp_chr <- as.integer(cnv_block[i,1 ] )
    hic_matrix_of_cnv_block[[i]] <- FindLocalHicMatrix(region1 = cnv_block[i, ],
                                                       region2 = cnv_block[i, ],
                                                       matrix = as.matrix( hic_matrix_list[[tmp_chr ]][[ tmp_chr ]])
                                                       )
  }
  if (methods == 'median' ) {
    hic_matrix_of_cnv_block_median <- sapply(hic_matrix_of_cnv_block, function (x) median( x[upper.tri(x, diag = T)], na.rm= T ) )
  }
  else if ( methods == 'mean') {
    hic_matrix_of_cnv_block_median <- sapply(hic_matrix_of_cnv_block, function (x) mean( x[upper.tri(x, diag = T)], na.rm= T ) )
  }
  
  return(  hic_matrix_of_cnv_block_median )
}

DiagHicScoreInCNVBlock <- function ( cnv_block, hic_matrix_list) {
  hic_matrix_of_cnv_block <- list()
  for ( i in 1:dim(cnv_block)[1]) {
    tmp_chr <- as.integer(cnv_block[i,1 ] )
    hic_matrix_of_cnv_block[[i]] <- FindLocalHicMatrix(region1 = cnv_block[i, ],
                                                       region2 = cnv_block[i, ],
                                                       matrix = as.matrix( hic_matrix_list[[tmp_chr ]][[ tmp_chr ]])
    )
  }
  hic_matrix_of_cnv_block_diag <- sapply(hic_matrix_of_cnv_block, function (x) median( diag(x), na.rm= T ) )
  return(  hic_matrix_of_cnv_block_diag )
}
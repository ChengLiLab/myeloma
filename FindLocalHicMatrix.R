# find local hic matrix
# 20160112
# wupz
FindLocalHicMatrix <- function(region1, region2, matrix , resolution = 200000) {
  # diag(matrix) <- diag(matrix)/2
  # find the start and end id
  # notice the boudary of chromosomes
  start1 <- floor( as.numeric(region1[2])/resolution ) + 1
  start2 <- floor( as.numeric(region2[2])/resolution ) + 1
  end1 <- ceiling( as.numeric(region1[3])/resolution ) 
  end2 <- ceiling( as.numeric(region2[3])/resolution ) 
  # return the matrix
  return(as.matrix(matrix[start1:end1, start2:end2] ) )
  
}

# test 
# FindLocalHicMatrix( RMPI_8226_cnv_block_200kb[1, ], RMPI_8226_cnv_block_200kb[1, ], raw_HiC_matrix[[1]][[1]][[1]], resolution = 2000000)
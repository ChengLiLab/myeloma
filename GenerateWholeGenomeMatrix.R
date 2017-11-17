#20151217
#wupz
GenerateWholeGenomeMatrix <- function( matrixList, cis = T ) {
  # Generate a whole genemo Hi-C matrix from a list of cis and trans matrix
  #
  # Args: 
  #   matrixList: a matrix list read from the disk
  # Returns:
  #   a whole genemo Hi-C matrix 
  dim_vector <- sapply( matrixList[[1]], dim)[2, ] 
  total_dim_size <- sum(dim_vector )  #get the total dim size
  wholeGenomeMatrix <- matrix(NA, nrow = total_dim_size, ncol = total_dim_size )
  #get the dim index of each chromosomes
  index_2 <- cumsum( dim_vector)
  index_1 <- c(1, index_2 + 1)
  for (i in 1:length(matrixList)) {
    for (j in 1:length(matrixList) ) {
      if( j < i) {
        wholeGenomeMatrix[index_1[i]:index_2[i], index_1[j]:index_2[j] ] <- t( as.matrix(matrixList[[j]][[i]] ) )
      }
      else if ( j == i) {
        if( cis == T) {
          wholeGenomeMatrix[index_1[i]:index_2[i], index_1[j]:index_2[j] ] <- as.matrix(matrixList[[i]][[i]] )
        }
        else if ( cis == F) {
          wholeGenomeMatrix[index_1[i]:index_2[i], index_1[j]:index_2[j] ] <- matrix(0, nrow = dim(as.matrix(matrixList[[i]][[i]] ))[1], 
                                                                                     ncol = dim(as.matrix(matrixList[[i]][[i]] ))[2])
        }
      }
      else {
        wholeGenomeMatrix[index_1[i]:index_2[i], index_1[j]:index_2[j] ] <- as.matrix(matrixList[[i]][[j]] )
      }
    }
  }
  return( list(dim_vector, wholeGenomeMatrix) )
  
}
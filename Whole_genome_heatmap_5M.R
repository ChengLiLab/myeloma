source("../scripts/GenerateWholeGenomeMatrix.R")
Whole_genome_heatmap <- list()
for (i in 1:4 ) {
  for ( j in 1:24 ) {
    tmp <- raw_HiC_matrix_5M[[i]][[j]][[j]]
    tmp_0_matrix <- matrix(data = 0, nrow = dim( raw_HiC_matrix_5M[[i]][[j]][[j]])[1] , ncol = dim( raw_HiC_matrix_5M[[i]][[j]][[j]])[2] )
    raw_HiC_matrix_5M[[i]][[j]][[j]] <- tmp_0_matrix
  }
  Whole_genome_heatmap[[i]] <- GenerateWholeGenomeMatrix(matrixList = raw_HiC_matrix_5M[[i]])
}

save(Whole_genome_heatmap, file = "../preprocessingData/Whole_genome_heatmap_5M.RData")
Whole_genome_heatmap[[1]][[1]]

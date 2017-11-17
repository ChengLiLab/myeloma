getLocalHeatmap <- function( query_sites, hic_matrix) {
  local_heatmap <- matrix(NA, nrow = 11, ncol = 11)
  chr1 <- query_sites[1]
  chr2 <- query_sites[2]
  bin1 <- query_sites[3]
  bin2 <- query_sites[4]
  chr1_chr2_heatmap <- hic_matrix[[chr1]][[chr2]]
  chr1_chr2_dims <- dim(chr1_chr2_heatmap)
  # get seq 1 in chr 1
  if(bin1 <= 5 ) {
    seq1 <- 1:(bin+5)
    if ( bin 2 <= 5 ) {
      seq2 <- 1:(bin+5) 
      tmp_local_heatmap <- chr1_chr2_heatmap[seq1, seq2]
      local_heatmap[(7 - bin1):11, (7 - bin2):11] <- tmp_local_heatmap
    }
    else if ( bin2 > (chr1_chr2_dims[2] -5) ) {
      seq2 <- (bin2-5):chr1_chr2_dims[2]
      tmp_local_heatmap <- chr1_chr2_heatmap[seq1, seq2]
      local_heatmap[(7 - bin1):11, 1:(chr1_chr2_dims[2] - bin2 + 4)] <- tmp_local_heatmap
    }
    else {
      seq2 <- (bin2-5):(bin2 + 5)
      tmp_local_heatmap <- chr1_chr2_heatmap[seq1, seq2]
      local_heatmap[(7 - bin1):11, ] <- tmp_local_heatmap
    }
  }
  # 
  else if ( bin1 > (chr1_chr2_dims[1] -5) ) {
    seq1 <- (bin1-5):chr1_chr2_dims[1]
  }
  else { seq1 <- (bin1 - 5):(bin1 + 5) }
  # get seq 2 in chr 2
  if( bin2 <= 5 ) {
    seq2 <- 1:(bin+5)
  }
  else if ( bin2 > (chr1_chr2_dims[2] -5) ) {
    seq2 <- (bin2-5):chr1_chr2_dims[2]
  }
  else { seq2 <- (bin2 - 5):(bin2 + 5) }
  tmp_local_heatmap <- chr1_chr2_heatmap[seq1, seq2]
  
  return(tmp_local_heatmap )
  
}

# test
getLocalHeatmap ( ctx_score[1, ], hic_matrix = raw_HiC_matrix[[1]])



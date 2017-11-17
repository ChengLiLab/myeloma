# input a pair of translocation sites and Hi-C matrix, get its hi-c interaction score
getHicScore <- function(query, hic_matrix) {
  chr1 <- as.integer(strsplit(query[1], split = "chr")[[1]][2])
  bin1 <- ceiling( as.numeric(query[2]) / 200000)
  chr2 <- as.integer(strsplit(query[5], split = "chr")[[1]][2])
  bin2 <- ceiling( as.numeric(query[6]) / 200000)
  if ( chr2 < chr1) {
    tmp  <- chr2
    chr2 <- chr1
    chr1 <- tmp
    tmp <- bin2
    bin2 <- bin1
    bin1 <- tmp
    
  }
  score_1 <- hic_matrix[[chr1]][[chr2]][bin1, bin2]
  #score_1 <- hic_matrix[[chr1]][[chr2]][bin1, bin2]
  #try(score_2 <- hic_matrix[[chr1]][[chr2]][bin1 +1 , bin2 +1], silent = TRUE)
  
  #if ( !exists("score_2") ) {
  # score_2 <- 0
  #}
  return(c(chr1, chr2, bin1, bin2, score_1))
}
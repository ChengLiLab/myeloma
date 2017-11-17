# combine two chromosome interaction matrix
# 20160106
# wupz
CombineTwoChromosomes <- function ( raw_matrix , chr_i, chr_j) {
  combine_chri_chrj <- cbind( raw_matrix[[chr_i]][[chr_i]],raw_matrix[[chr_i]][[chr_j]] )
  combine_chri_chrj <- rbind( combine_chri_chrj, cbind(t(raw_matrix[[chr_i]][[chr_j]]), raw_matrix[[chr_j]][[chr_j]] ) )
  print( dim(combine_chri_chrj) )
  return( combine_chri_chrj)
}



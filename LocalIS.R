# find local insulation score
# wupz
# 20160128
LocalIS <- function( bed_vec, insulation_score, flank = 25) {
  # get the resolution of insulation_score
  resolution <- insulation_score[1, 'end'] - insulation_score[1, 'start']
  # add 0 data to the head and end of the table to avoid border problem
  tmp_0 <- matrix(NA, nrow = flank, ncol = dim(insulation_score)[2] )
  tmp_ins_matrix <- rbind(tmp_0, as.matrix(insulation_score), tmp_0)
  # find the region id of input region in insulation_score +- flank
  start_neibour_id <- floor(as.numeric(bed_vec[2]) / resolution) + 1 + flank
  #end_neibour_id <- floor(as.numeric(bed_vec[3]) / resolution) + 1 + flank
  start_neibour_insulation_score <- tmp_ins_matrix[(start_neibour_id - flank):(start_neibour_id + flank), 'insulationScore']
  #end_neibour_insulation_score <- tmp_ins_matrix[(end_neibour_id - flank):(end_neibour_id + flank), 'insulationScore']
  return( start_neibour_insulation_score )
}



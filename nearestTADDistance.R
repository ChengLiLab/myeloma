# the nearest distance of a sites from TAD boundary
nearestTADDistance <- function(query_Sites, subject_Sites) {
  query <- IRanges( as.numeric( query_Sites[, 2]) , as.numeric(query_Sites[, 2]) )
  subject <- IRanges( subject_Sites[, 2] -1 , subject_Sites[,2] -1 )
  nearest_distance <- distanceToNearest(query, subject )@elementMetadata[[1]]
  return(nearest_distance)
}
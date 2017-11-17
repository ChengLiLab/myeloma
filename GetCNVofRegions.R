## input a region and a cnv file
## retrive the cnv of this region
## wupz
## 20170510
GetCNVofRegions <- function (cnv_files, region_bed, methods = "median") {
  if( !is.numeric(i) ) {
    tmp_chr <- as.integer(strsplit(region_bed[1], split = "chr")[[1]][[2]])
  }
  else {
    tmp_chr <- region_bed[1]
  }
  tmp_resolution <- as.numeric(cnv_files[2, 2]) - as.numeric(cnv_files[2, 1])
  tmp_region_start_id <- ceiling(as.numeric(region_bed[2])/tmp_resolution)
  tmp_region_end_id <- floor(as.numeric(region_bed[3])/tmp_resolution)
  tmp_region_id <- tmp_region_start_id:tmp_region_end_id
  if ( methods == "median") {
    regions_cnv <- floor( median(cnv_files[cnv_files[, 1] == tmp_chr, ][tmp_region_id, 5] ) )
  } else if ( methods == "mean" ) {
    regions_cnv <- mean(cnv_files[cnv_files[, 1] == tmp_chr, ][tmp_region_id, 5] )
  }
  return(regions_cnv)
} 

# test
# GetCNVofRegions(cnv_files = U266_CNV_40kb, region_bed = U266_TAD_filter_bed[1, ] )

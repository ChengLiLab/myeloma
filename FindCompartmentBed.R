FindCompartment <- function( hic_pca )  {
  # if ratio equals to -1, set the copy number to 0
  resolution = hic_pca$start[2] - hic_pca$start[1]
  hic_pca$compartment[is.na(hic_pca$compartment)] <- 0
  j = 1
  while( hic_pca[j, 4] == 0 ) {
    j <- j + 1
  }
  tmp_compartment_id <- 1
  # save the results in a bed format object
  compartment_bed <- data.frame(chr = i , 
                                start = (j-1) * resolution, # 0-base
                                end = j*resolution,
                                name = paste("Compartment", tmp_compartment_id, sep = "_"),
                                score = hic_pca[j, 5],
                                strand = '.', 
                                stringsAsFactors = FALSE
  )
  # print(j)
  for( k in (j+1):dim(hic_pca)[1] ) {
    if (hic_pca[k, 5] == 0) {
      print(k)
    }
    else if( hic_pca[k, 1] != hic_pca[k-1, 1] | hic_pca[k, 5] != hic_pca[k-1, 5] ) {
      tmp_compartment_id <- tmp_compartment_id + 1
      print( tmp_compartment_id )
      add_compartment <- data.frame(chr =  hic_pca[k, 1] , 
                                    start = as.numeric(hic_pca[k, 2] ), # change from 1 based to 0-based
                                    end = as.numeric(hic_pca[k, 3] ),
                                    name = paste('Compartment',tmp_compartment_id , sep = '_'),
                                    score = hic_pca[k, 5],
                                    strand = '.', 
                                    stringsAsFactors = FALSE
                                    )
      compartment_bed <- rbind( compartment_bed, add_compartment)
    }
    else if ( hic_pca[i, 5] == hic_pca[i-1, 5] ) {
      compartment_bed[dim(compartment_bed)[1], 3] <- hic_pca[k, 3]
    }
  }
  return(compartment_bed )
}

# test
test <- FindCompartment(U266_compartment_MobI)

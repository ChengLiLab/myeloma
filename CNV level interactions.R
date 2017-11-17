# CNV level interactions
# 20170411
# wupz

# read data
U266_CNV_40kb <- read.table("/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_40kb/U266-merged.sorted.bam.dedup.bam_ratio.txt",
                            stringsAsFactors = F, header = T)
head(U266_CNV_40kb)
dim(U266_CNV_40kb)
U266_CNV_40kb_block <- FindCNVBlock(freec_ratio = U266_CNV_40kb)
head(U266_CNV_40kb_block)
U266_CNV_40kb_block[U266_CNV_40kb_block[, 1] == 5,]

# calculate CNV level interaction matrix
tmp_dim <- dim(U266_CNV_40kb_block[U266_CNV_40kb_block[, 1] != "Y" & U266_CNV_40kb_block[, 1] != "X", ])[1] # ignore chrX and chrY
U266_CNV_40kb_block_hic_matrix <- matrix(0, nrow = tmp_dim, ncol = tmp_dim )
tmp_chr_length <- 0
for (tmp_chr in paste0("chr", 1:22) ) {
  tmp_u266_hic <- read.table( paste( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/ice_normalization/", 
                                     tmp_chr, "_40k_normalized_matrix.txt", sep = "") )
  tmp_u266_hic <- as.matrix(tmp_u266_hic)
  diag(tmp_u266_hic) <- diag(tmp_u266_hic) *2
  tmp_cnv_file <- U266_CNV_40kb_block[U266_CNV_40kb_block[, 1] == strsplit(tmp_chr, "chr")[[1]][[2]], ]
  
  for(i in 1:dim(tmp_cnv_file)[1] ) {
    for (j in i:dim(tmp_cnv_file)[1] ) {
      U266_CNV_40kb_block_hic_matrix[ tmp_chr_length + i, tmp_chr_length + j] <- median( as.vector(as.numeric( FindLocalHicMatrix(region1 = tmp_cnv_file[i, ], 
                                                                                                region2 = tmp_cnv_file[j, ], 
                                                                                                matrix = tmp_u266_hic, resolution = 40000) ) ), 
                                                       na.rm = T
      )
    }
  }
  tmp_chr_length <- tmp_chr_length + dim(tmp_cnv_file)[1]
}

# stat for the gain and loss interactions
# found TAD interactions seperate by another TAD
dim(U266_CNV_40kb_block_hic_matrix)[1] # [1] 401
tmp_U266_CNV_40kb_block_hic_matrix <- U266_CNV_40kb_block_hic_matrix[ -c(dim(U266_CNV_40kb_block_hic_matrix)[1], dim(U266_CNV_40kb_block_hic_matrix)[1]-1), -c(1,2)]
dim(tmp_U266_CNV_40kb_block_hic_matrix)
tmp_U266_CNV_40kb_block_hic_matrix_vector <- diag(tmp_U266_CNV_40kb_block_hic_matrix)
# plot interaction between TADs steps by 2
plot(density(tmp_U266_CNV_40kb_block_hic_matrix_vector))
# calculate distance between every 2 CNVs
tmp_U266_CNV_median_point <- U266_CNV_40kb_block[, 2]/2 +  U266_CNV_40kb_block[, 3]/2
tmp_U266_CNV_median_point <- tmp_U266_CNV_median_point[1:tmp_dim]
length(tmp_U266_CNV_median_point)
tmp_U266_CNV_40kb_block_distance <- tmp_U266_CNV_median_point[-c(1,2)] -  tmp_U266_CNV_median_point[-c(length(tmp_U266_CNV_median_point)[1], length(tmp_U266_CNV_median_point)[1] -1)]
length(tmp_U266_CNV_40kb_block_distance)
tmp_U266_TAD_filter_distance_vector <- diag(tmp_U266_TAD_filter_distance_matrix)
plot(x = tmp_U266_CNV_40kb_block_distance, y = tmp_U266_CNV_40kb_block_hic_matrix_vector, xlim = c(0, 100000000) )

png(filename = "U266_chr1_TAD_lever_interactions_vs_distance.png", width = 1024, height = 1024)
tmp_cnv <- U266_CNV_40kb_block[1:tmp_dim, 5]
cnv_gain_loss <- numeric() # used for store the gain or loss state
for (i in 1:(length(tmp_cnv)-2) ) {
  if ( (tmp_cnv[i+1] < tmp_cnv[i]) & (tmp_cnv[i+1] < tmp_cnv[i+2]) ) {
    cnv_gain_loss[i-1] <- "Loss"
  } else if ( (tmp_cnv[i+1] > tmp_cnv[i]) & (tmp_cnv[i+1] > tmp_cnv[i+2]) ) {
    cnv_gain_loss[i] <- "Gain"
  } else {
    cnv_gain_loss[i] <- "Others"
  }
}


cnv_gain_loss_col <- ifelse(cnv_gain_loss == "Gain", yes = "red", no = ifelse(cnv_gain_loss == "Loss", "blue", "green") )
plot(x=tmp_U266_CNV_40kb_block_distance, y = tmp_U266_CNV_40kb_block_hic_matrix_vector, col = cnv_gain_loss_col, xlim = c(0, 100000000), 
     xlab = "Distance", ylab = "Interactions")
boxplot(tmp_U266_CNV_40kb_block_hic_matrix_vector[tmp_U266_CNV_40kb_block_hic_matrix_vector!=0]~cnv_gain_loss[tmp_U266_CNV_40kb_block_hic_matrix_vector!=0], ylim = c(0, 50))
dev.off()




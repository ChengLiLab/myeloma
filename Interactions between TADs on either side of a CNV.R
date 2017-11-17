# Interactions between TADs on either side of a CNV
# 201170412
# wupz

# 1.1 read CNV data
U266_CNV_40kb <- read.table("/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_40kb/U266-merged.sorted.bam.dedup.bam_ratio.txt",
                            stringsAsFactors = F, header = T)
head(U266_CNV_40kb)
dim(U266_CNV_40kb)
U266_CNV_40kb_block <- FindCNVBlock(freec_ratio = U266_CNV_40kb)
head(U266_CNV_40kb_block )
dim(U266_CNV_40kb_block )

# 1.2 find the gain or loss cnv regions
tmp_dim <- dim(U266_CNV_40kb_block)[1] # 
tmp_cnv <- U266_CNV_40kb_block$score
length(tmp_cnv)
cnv_gain_loss_state <- vector() # used for store the state of gain or loss

for (i in 1:(length(tmp_cnv)-2) ) {
  if ( (U266_CNV_40kb_block$chrom[i+1] ==  U266_CNV_40kb_block$chrom[i]) & (U266_CNV_40kb_block$chrom[i+1] ==  U266_CNV_40kb_block$chrom[i+2]) ) { # the same chrom
    if ( (tmp_cnv[i+1] < tmp_cnv[i]) & (tmp_cnv[i+1] < tmp_cnv[i+2]) ) { # loss fragment
      cnv_gain_loss_state[i] <- "Loss"
    } else if ( (tmp_cnv[i+1] > tmp_cnv[i]) & (tmp_cnv[i+1] > tmp_cnv[i+2]) ) { # gain fragment
      cnv_gain_loss_state[i] <- "Gain"
    } else {
      cnv_gain_loss_state[i] <- "Others"
    }
  } else {
    cnv_gain_loss_state[i] <- "Others" #  start or end of a chromo
  }
}
length(cnv_gain_loss_state) # [1] 428
cnv_gain_loss_state <- c("Others", cnv_gain_loss_state, "Others")
U266_CNV_40kb_block <- cbind(U266_CNV_40kb_block, cnv_gain_loss_state = cnv_gain_loss_state)
class(U266_CNV_40kb_block)
dim(U266_CNV_40kb_block)
U266_CNV_40kb_block$chrom[U266_CNV_40kb_block$chrom=="X"] <- 23
U266_CNV_40kb_block$chrom[U266_CNV_40kb_block$chrom=="Y"] <- 24
sum(U266_CNV_40kb_block$chromEnd - U266_CNV_40kb_block$chromStart <= 120000)
png(filename = "Length of U266 CNV blocks.png", width = 1024, height = 1024)
par(cex = 3, lwd = 3, mgp = c(2.5,1,0))
plot(density(U266_CNV_40kb_block$chromEnd-U266_CNV_40kb_block$chromStart), log = "x", 
     xlab = "Length of CNV blocks (bp)", main = "Length of U266 CNV blocks")
# lines(density(RPMI8226_CNV_40kb_block$chromEnd-RPMI8226_CNV_40kb_block$chromStart), col = "red")
dev.off()
write.table(x = U266_CNV_40kb_block, file = "U266_CNV_40kb_block.bed", quote = F, sep = "\t", row.names = F, col.names = F)

# 1.3 For each TAD, calculate region starts and ends and CNVs of this TAD
## use 1...23 as chromosome names
tmp_U266_CNV_40kb <- U266_CNV_40kb
tmp_U266_CNV_40kb$Chromosome[U266_CNV_40kb$Chromosome == "X"] <- 23
tmp_U266_CNV_40kb$Chromosome[U266_CNV_40kb$Chromosome == "Y"] <- 24
# read insulation score data
U266_insulation_score <- list()
U266_insulation_score_tad <- list()
for ( i in 1:23 ) {
  tmp_chr <- paste0("chr", i)
  U266_insulation_score[[i]] <- read.table(paste( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/TAD_boundary/", 
                                                 tmp_chr, "_", tmp_chr, "_40k_normalmatrix.txt.is1000001.ids240001.insulation", sep = ""), header = T)
  U266_insulation_score_tad[[i]] <- read.table(paste( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/TAD_boundary/", 
                                                     tmp_chr, "_", tmp_chr, "_40k_normalmatrix.txt.is1000001.ids240001.insulation.boundaries", sep = ""), header = T)
}
U266_TAD_filter_bed <- vector()
U266_TAD_bed <- vector()
for ( i in 1:23  ) {
  tmp_chr <- paste0("chr", i)
  tmp_u266_insulation_score <- U266_insulation_score[[i]] 
  tmp_u266_insulation_score_tad <- U266_insulation_score_tad[[i]]
  # construct TAD file with bed format
  tmp_U266_TAD_bed <- cbind(chr = i, 
                        TAD_start = tmp_u266_insulation_score_tad$end[-length(tmp_u266_insulation_score_tad$end)] - 1, 
                        TAD_end = tmp_u266_insulation_score_tad$end[-1] - 40001)
  dim(tmp_U266_TAD_bed) # Check numbers
  tmp_U266_TAD_bed <- cbind(tmp_U266_TAD_bed, apply(X = tmp_U266_TAD_bed, MARGIN = 1, FUN = GetCNVofRegions, cnv_files = tmp_U266_CNV_40kb))
  U266_TAD_bed <- rbind(U266_TAD_bed, tmp_U266_TAD_bed)
  # keep TAD which length over 160 kb
  tmp_U266_TAD_filter_bed <- tmp_U266_TAD_bed[ as.numeric(tmp_U266_TAD_bed[, 3]) - as.numeric(tmp_U266_TAD_bed[, 2]) >= 160000, ]
  dim(tmp_U266_TAD_filter_bed) # Check numbers
  # calculate CNV at TAD level
  head(tmp_U266_TAD_filter_bed )
  U266_TAD_filter_bed <- rbind(U266_TAD_filter_bed, tmp_U266_TAD_filter_bed)
}
# check data
dim( U266_TAD_bed )
dim( U266_TAD_filter_bed )
png(filename = "U266 Length of raw TADs.png", width = 1024, height = 1024)
par(cex = 3, lwd = 3, mgp = c(2.5,1,0))
plot(table(as.numeric(U266_TAD_bed[, 3])-as.numeric(U266_TAD_bed[, 2])), xlim = c(0, 3000000), type = "l",
     xlab = "Length of \"raw\" TADs (bp)", main = "U266, distribution of TADs' length", ylab = "Counts")
dev.off()

png(filename = "U266 Length of filtered TADs.png", width = 1024, height = 1024)
par(cex = 3, lwd = 3, mgp = c(2.5,1,0))
plot(table(as.numeric(U266_TAD_filter_bed[, 3])-as.numeric(U266_TAD_filter_bed[, 2])), xlim = c(0, 3000000), type = "l",
     xlab = "Length of \"filtered\" TADs (bp)", main = "U266, distribution of TADs' length", ylab = "Counts")
dev.off()
write.table( x = U266_TAD_bed, file = "U266_TAD_raw.bed", quote = F, sep = "\t", col.names = F, row.names = F )
write.table( x = U266_TAD_filter_bed, file = "U266_TAD_filter.bed", quote = F, sep = "\t", col.names = F, row.names = F )

## found local intra-chromosome interactions in each TAD
# read Hi-C data of U266, 40kb
U266_hic_ice_matrix <- list()
U266_hic_raw_matrix <- list()
for (i in 1:23) {
  print(paste0("Read raw U266 Hi-C matrix of chr ", i) )
  U266_hic_ice_matrix[[i]] <- list()
  U266_hic_ice_matrix[[i]][[i]] <- read.table(paste0( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/ice_normalization/chr", 
                                                  i, "_40k_normalized_matrix.txt") )
  print(paste0("Read raw U266 Hi-C matrix of chr ", i) )
  U266_hic_raw_matrix[[i]] <- list()
  U266_hic_raw_matrix[[i]][[i]] <- read.table(paste0( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/raw/chr", 
                                                i, "_40k_rawmatrix.txt") )
}
# calculate TAD-TAD level interaction, use median interaction as represent
U266_TAD_filter_HiC_matrix_raw <- list()
U266_TAD_filter_HiC_matrix_ice <- list()
for (i in 1:23) {
  tmp_U266_TAD_filter_bed <- U266_TAD_filter_bed[U266_TAD_filter_bed[, 1] == i, ]
  U266_TAD_filter_HiC_matrix_raw[[i]] <- matrix(0, nrow = dim(tmp_U266_TAD_filter_bed)[1],  ncol = dim(tmp_U266_TAD_filter_bed)[1])
  U266_TAD_filter_HiC_matrix_ice[[i]] <- matrix(0, nrow = dim(tmp_U266_TAD_filter_bed)[1],  ncol = dim(tmp_U266_TAD_filter_bed)[1])
  tmp_u266_hic_raw <- U266_hic_raw_matrix[[i]][[i]]
  tmp_u266_hic_ice <- U266_hic_ice_matrix[[i]][[i]]
  diag(tmp_u266_hic_raw) <- diag(tmp_u266_hic_raw)/2
  diag(tmp_u266_hic_ice) <- diag(tmp_u266_hic_ice)/2
  print("Calculate TAD-TAD interaction")
  for (j in 1:dim(tmp_U266_TAD_filter_bed)[1] ) {
    for (k in i:dim(tmp_U266_TAD_filter_bed)[1] ) {
      U266_TAD_filter_HiC_matrix_raw[[i]][j, k] <- median( as.vector(as.numeric( FindLocalHicMatrix(region1 = tmp_U266_TAD_filter_bed[j, ], 
                                                                                           region2 = tmp_U266_TAD_filter_bed[k, ], 
                                                                                           matrix = tmp_u266_hic_raw, 
                                                                                           resolution = 40000) ) ), 
                                                  na.rm = T)
      U266_TAD_filter_HiC_matrix_ice[[i]][j, k] <- median( as.vector(as.numeric( FindLocalHicMatrix(region1 = tmp_U266_TAD_filter_bed[j, ], 
                                                                                                    region2 = tmp_U266_TAD_filter_bed[k, ], 
                                                                                                    matrix = tmp_u266_hic_ice, 
                                                                                                    resolution = 40000) ) ), 
                                                           na.rm = T)
    }
  }
  
}

U266_local_intra_hic_matrix_raw <- 


# 1.4 U266, Read up- and down-stream TAD regions and calculate its interactions
# 1.4.1 Read files of the upstream or downstream closest TADs to the CNVs
U266_CNV_40kb_upTAD <- read.table("/lustre/user/liclab/wupz/dosageEffect/201703_nc_revise/FISH_TAD_CNV_interactions/U266_CNV_nearestUpstream_TAD.bed", 
                                  stringsAsFactors = F)
U266_CNV_40kb_downTAD <- read.table("/lustre/user/liclab/wupz/dosageEffect/201703_nc_revise/FISH_TAD_CNV_interactions/U266_CNV_nearestDownstream_TAD.bed", 
                                  stringsAsFactors = F)
U266_CNV_40kb_nearest_TAD <- cbind(U266_CNV_40kb_upTAD, U266_CNV_40kb_downTAD[, 8:12])
head(U266_CNV_40kb_nearest_TAD)
dim(U266_CNV_40kb_nearest_TAD)
# 1.4.2 calculate the upstreamTAD-downstreamTAD interactions
U266_CNV_40kb_nearest_TAD_interactions <- numeric(length = dim(U266_CNV_40kb_nearest_TAD)[1])
tmp_chr <- "chr1"
tmp_u266_hic <- read.table( paste( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/ice_normalization/", 
                                   tmp_chr, "_40k_normalized_matrix.txt", sep = "") )
for( i in 1:dim(U266_CNV_40kb_nearest_TAD)[1] ) {
  # change hi-c data according to chromosomes
  if ( tmp_chr != paste0("chr", U266_CNV_40kb_upTAD[i, 1] ) ) {
    tmp_chr <- paste0("chr", U266_CNV_40kb_upTAD[i, 1] )
    tmp_u266_hic <- read.table( paste( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/ice_normalization/", 
                                       tmp_chr, "_40k_normalized_matrix.txt", sep = "") )
  }
  if( U266_CNV_40kb_upTAD[i, 8]!= "." &  U266_CNV_40kb_downTAD[i, 8] != "." ) { # not include the border
  U266_CNV_40kb_nearest_TAD_interactions[i] <- median( as.vector( as.numeric( FindLocalHicMatrix(region1 = U266_CNV_40kb_upTAD[i, 8:10], 
                                                                                                 region2 = U266_CNV_40kb_downTAD[i, 8:10], 
                                                                                                 matrix = tmp_u266_hic, resolution = 40000) ) ) )
  }
}
## 1.4.3 Generate two random TADs of the same CNV block and find its interactions
### 1.4.3.1 choose large CNV blocks
tmp_index <- U266_CNV_40kb_block$chromEnd - U266_CNV_40kb_block$chromStart >= 5000000
sum( tmp_index ) # [1] 121
head( U266_CNV_40kb_block[tmp_index,] )
### 1.4.3.2 select two TADs from the chose CNV blocks
for( tmp_chr in U266_CNV_40kb_block[tmp_index, 1] ) {
  
}
## 1.4.3 color calculation and distance calculation
tmp_col <- ifelse(U266_CNV_40kb_upTAD[, 7] == "Gain", "red", ifelse(U266_CNV_40kb_upTAD[, 7] == "Loss", "blue", "green"))
tmp_up_distance <- abs(U266_CNV_40kb_upTAD[, 12])
tmp_down_distance <- abs(U266_CNV_40kb_downTAD[, 12]) # tmp_distance < 260000 & tmp_distance_2 < 260000 & 
tmp_index <- U266_CNV_40kb_downTAD[, 8] != "." & U266_CNV_40kb_upTAD[, 8] != "."
tmp_distance_3 <- apply(U266_CNV_40kb_downTAD[tmp_index, c(9,10)], 1, mean) - apply(U266_CNV_40kb_upTAD[tmp_index, c(9,10)], 1, mean)
## 1.4.4 Statistic figures 
### Interactions vs. distance
png(filename = "U266_Interactions_upTAD_downTAD_of_CNV_vs_distance.png", width = 1024, height = 1024)
par(cex = 2.5, mgp = c(1, 0.5, 0))
plot(U266_CNV_40kb_nearest_TAD_interactions[tmp_index] ~ tmp_distance_3, col = tmp_col[tmp_index])
dev.off()
### Interactions (all) vs. cnv state
png(filename = "U266_Interactions_upTAD_downTAD_of_CNV_vs_cnvState.png", width = 1024, height = 1024)
par(cex = 2.5, mgp = c(1, 0.5, 0))
boxplot( U266_CNV_40kb_nearest_TAD_interactions[tmp_index] ~ U266_CNV_40kb_upTAD[tmp_index, 7],
         col = c("red", "blue", "green"))
dev.off()
### Interactions (non-zero) vs. cnv state
png(filename = "U266_nonZero_Interactions_upTAD_downTAD_of_CNV_vs_cnvState.png", width = 1024, height = 1024)
par(cex = 2.5, mgp = c(1, 0.5, 0))
boxplot( U266_CNV_40kb_nearest_TAD_interactions[tmp_index][U266_CNV_40kb_nearest_TAD_interactions[tmp_index]!=0] ~ U266_CNV_40kb_upTAD[tmp_index, 7][U266_CNV_40kb_nearest_TAD_interactions[tmp_index]!=0],
         col = c("red", "blue", "green"))
dev.off()
tmp_interactions <- U266_CNV_40kb_nearest_TAD_interactions[tmp_index][U266_CNV_40kb_nearest_TAD_interactions[tmp_index]!=0]
tmp_cnv_state <- U266_CNV_40kb_upTAD[tmp_index, 7][U266_CNV_40kb_nearest_TAD_interactions[tmp_index]!=0]
t.test(tmp_interactions[tmp_cnv_state=="Loss"], tmp_interactions[tmp_cnv_state=="Others"])
t.test(tmp_interactions[tmp_cnv_state=="Gain"], tmp_interactions[tmp_cnv_state=="Others"])
### Distance vs. cnv state
png(filename = "U266_nonZero_Interactions_upTAD_downTAD_of_CNV_vs_cnvState.png", width = 1024, height = 1024)
par(cex = 2.5, mgp = c(1, 0.5, 0))
boxplot(tmp_distance_3 ~ U266_CNV_40kb_upTAD[tmp_index, 7], col = c("red", "blue", "green"))
dev.off()





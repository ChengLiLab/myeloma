# 2 RPMI8226, Interactions between TADs on either side of a CNV
# 2.1 read CNV data
RPMI8226_CNV_40kb <- read.table("/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_40kb/8226-merged.sorted.bam.dedup.bam_ratio.txt",
                                stringsAsFactors = F, header = T)
head(RPMI8226_CNV_40kb)
dim(RPMI8226_CNV_40kb)
RPMI8226_CNV_40kb_block <- FindCNVBlock(freec_ratio = RPMI8226_CNV_40kb)
dim(RPMI8226_CNV_40kb_block) # 200kb: [1] 272   6, 40kb: [1] 596   6
RPMI8226_CNV_40kb_block[RPMI8226_CNV_40kb_block[, 1] == 5,]

# 2.2 find the gain or loss cnv regions
tmp_dim <- dim(RPMI8226_CNV_40kb_block)[1] # 
tmp_cnv <- RPMI8226_CNV_40kb_block$score
length(tmp_cnv)
cnv_gain_loss_state <- vector() # used for store the state of gain or loss

for (i in 1:(length(tmp_cnv)-2) ) {
  if ( (RPMI8226_CNV_40kb_block$chrom[i+1] ==  RPMI8226_CNV_40kb_block$chrom[i]) & (RPMI8226_CNV_40kb_block$chrom[i+1] ==  RPMI8226_CNV_40kb_block$chrom[i+2]) ) {
    if ( (tmp_cnv[i+1] < tmp_cnv[i]) & (tmp_cnv[i+1] < tmp_cnv[i+2]) ) {
      cnv_gain_loss_state[i] <- "Loss"
    } else if ( (tmp_cnv[i+1] > tmp_cnv[i]) & (tmp_cnv[i+1] > tmp_cnv[i+2]) ) {
      cnv_gain_loss_state[i] <- "Gain"
    } else {
      cnv_gain_loss_state[i] <- "Others"
    }
  } else {
    cnv_gain_loss_state[i] <- "Others"
  }
}
length(cnv_gain_loss_state) # [1] 594
cnv_gain_loss_state <- c("Others", cnv_gain_loss_state, "Others")
RPMI8226_CNV_40kb_block <- cbind(RPMI8226_CNV_40kb_block, cnv_gain_loss_state = cnv_gain_loss_state)
class( RPMI8226_CNV_40kb_block )
dim( RPMI8226_CNV_40kb_block )
RPMI8226_CNV_40kb_block$chrom[RPMI8226_CNV_40kb_block$chrom=="X"] <- 23
RPMI8226_CNV_40kb_block$chrom[RPMI8226_CNV_40kb_block$chrom=="Y"] <- 24
write.table(x = RPMI8226_CNV_40kb_block, file = "RPMI8226_CNV_40kb_block.bed", quote = F, row.names = F, col.names = F)
# 2.3 For each TAD, calculate region starts and ends and CNVs
## use 1...23 as chromosome names
tmp_RPMI8226_CNV_40kb <- RPMI8226_CNV_40kb
tmp_RPMI8226_CNV_40kb$Chromosome[RPMI8226_CNV_40kb$Chromosome == "X"] <- 23
tmp_RPMI8226_CNV_40kb$Chromosome[RPMI8226_CNV_40kb$Chromosome == "Y"] <- 24
# 
RPMI8226_TAD_filter_bed <- vector()
for ( tmp_chr in paste0("chr", 1:23) ) {
  tmp_RPMI8226_insulation_score <- read.table(paste( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release10/resolution_40k/cis/TAD_boundary/", 
                                                     tmp_chr, "_", tmp_chr, "_40k_normalmatrix.txt.is1000001.ids240001.insulation", sep = ""), header = T)
  tmp_RPMI8226_insulation_score_tad <- read.table(paste( "/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release10/resolution_40k/cis/TAD_boundary/", 
                                                         tmp_chr, "_", tmp_chr, "_40k_normalmatrix.txt.is1000001.ids240001.insulation.boundaries", sep = ""), header = T)
  head(tmp_RPMI8226_insulation_score)
  head(tmp_RPMI8226_insulation_score_tad)
  # construct TAD file with bed format
  tmp_RPMI8226_TAD_bed <- cbind(chr = tmp_chr, 
                                TAD_start = tmp_RPMI8226_insulation_score_tad$end[-length(tmp_RPMI8226_insulation_score_tad$end)] - 1, 
                                TAD_end = tmp_RPMI8226_insulation_score_tad$end[-1] - 40001)
  dim(tmp_RPMI8226_TAD_bed) # Check numbers
  # keep TAD which length over 160 kb
  tmp_RPMI8226_TAD_filter_bed <- tmp_RPMI8226_TAD_bed[ as.numeric(tmp_RPMI8226_TAD_bed[, 3]) - as.numeric(tmp_RPMI8226_TAD_bed[, 2]) >= 160000, ]
  dim(tmp_RPMI8226_TAD_filter_bed) # Check numbers
  # calculate CNV at TAD level
  tmp_RPMI8226_TAD_filter_bed <- cbind(tmp_RPMI8226_TAD_filter_bed, apply(X = tmp_RPMI8226_TAD_filter_bed, MARGIN = 1, FUN = GetCNVofRegions, cnv_files = tmp_RPMI8226_CNV_40kb))
  head(tmp_RPMI8226_TAD_filter_bed )
  tmp_RPMI8226_TAD_filter_bed[, 1] <- strsplit(tmp_chr, "chr")[[1]][[2]]
  RPMI8226_TAD_filter_bed <- rbind(RPMI8226_TAD_filter_bed, tmp_RPMI8226_TAD_filter_bed)
  
}
## check data 
head(RPMI8226_TAD_filter_bed) 
tail(RPMI8226_TAD_filter_bed)
write.table(x = RPMI8226_TAD_filter_bed, file = "RPMI8226_TAD_filter.bed", quote = F, row.names = F)



# 3 GM12878, Interactions between TADs on either side of a CNV
# 3.1. read CNV data
GM12878_CNV_40kb <- read.table("/lustre/user/liclab/wupz/DNA_seq/CNV_final/GM12878_CNV/40kb/ERR091571.bam_ratio.txt",
                               stringsAsFactors = F, header = T)
head(GM12878_CNV_40kb)
dim(GM12878_CNV_40kb)
GM12878_CNV_40kb_block <- FindCNVBlock(freec_ratio = GM12878_CNV_40kb)
dim(GM12878_CNV_40kb_block) # [1] 238   6
GM12878_CNV_40kb_block[GM12878_CNV_40kb_block[, 1] == 5, ]

# find the gain or loss cnv regions
tmp_dim <- dim(GM12878_CNV_40kb_block)[1] # [1] 527 
tmp_cnv <- GM12878_CNV_40kb_block$score
length(tmp_cnv)
cnv_gain_loss_state <- vector() # used for store the state of gain or loss

for (i in 1:(length(tmp_cnv)-2) ) {
  if ( (GM12878_CNV_40kb_block$chrom[i+1] ==  GM12878_CNV_40kb_block$chrom[i]) & (GM12878_CNV_40kb_block$chrom[i+1] ==  GM12878_CNV_40kb_block$chrom[i+2]) ) {
    if ( (tmp_cnv[i+1] < tmp_cnv[i]) & (tmp_cnv[i+1] < tmp_cnv[i+2]) ) {
      cnv_gain_loss_state[i] <- "Loss"
    } else if ( (tmp_cnv[i+1] > tmp_cnv[i]) & (tmp_cnv[i+1] > tmp_cnv[i+2]) ) {
      cnv_gain_loss_state[i] <- "Gain"
    } else {
      cnv_gain_loss_state[i] <- "Others"
    }
  } else {
    cnv_gain_loss_state[i] <- "Others"
  }
}
length(cnv_gain_loss_state) # [1] 525
cnv_gain_loss_state <- c("Others", cnv_gain_loss_state, "Others")
GM12878_CNV_40kb_block <- cbind(GM12878_CNV_40kb_block, cnv_gain_loss_state = cnv_gain_loss_state)
class(GM12878_CNV_40kb_block)
dim(GM12878_CNV_40kb_block)
GM12878_CNV_40kb_block$chrom[GM12878_CNV_40kb_block$chrom=="X"] <- 23
GM12878_CNV_40kb_block$chrom[GM12878_CNV_40kb_block$chrom=="Y"] <- 24
write.table(x = GM12878_CNV_40kb_block, file = "GM12878_CNV_40kb_block.bed", quote = F, row.names = F, col.names = F)
# 3.2. For each TAD, calculate region starts and ends and CNVs
## use 1...23 as chromosome names
tmp_GM12878_CNV_40kb <- GM12878_CNV_40kb
tmp_GM12878_CNV_40kb$Chromosome[GM12878_CNV_40kb$Chromosome == "X"] <- 23
tmp_GM12878_CNV_40kb$Chromosome[GM12878_CNV_40kb$Chromosome == "Y"] <- 24
# 
GM12878_TAD_filter_bed <- vector()
for ( tmp_chr in paste0("chr", 1:23) ) {
  tmp_GM12878_insulation_score <- read.table(paste( "/lustre/user/liclab/lirf/Project/hic/2014/gm12878_hind3_dCTP_296mb/resolution_40k/cis/TAD_boundary/", 
                                                    tmp_chr, "_", tmp_chr, "_40k_normalmatrix.txt.is1000001.ids240001.insulation", sep = ""), header = T)
  tmp_GM12878_insulation_score_tad <- read.table(paste( "/lustre/user/liclab/lirf/Project/hic/2014/gm12878_hind3_dCTP_296mb/resolution_40k/cis/TAD_boundary/", 
                                                        tmp_chr, "_", tmp_chr, "_40k_normalmatrix.txt.is1000001.ids240001.insulation.boundaries", sep = ""), header = T)
  head(tmp_GM12878_insulation_score)
  head(tmp_GM12878_insulation_score_tad)
  # construct TAD file with bed format
  tmp_GM12878_TAD_bed <- cbind(chr = tmp_chr, 
                               TAD_start = tmp_GM12878_insulation_score_tad$end[-length(tmp_GM12878_insulation_score_tad$end)] - 1, 
                               TAD_end = tmp_GM12878_insulation_score_tad$end[-1] - 40001)
  dim(tmp_GM12878_TAD_bed) # Check numbers
  # keep TAD which length over 160 kb
  tmp_GM12878_TAD_filter_bed <- tmp_GM12878_TAD_bed[ as.numeric(tmp_GM12878_TAD_bed[, 3]) - as.numeric(tmp_GM12878_TAD_bed[, 2]) >= 160000, ]
  dim(tmp_GM12878_TAD_filter_bed) # Check numbers
  # calculate CNV at TAD level
  tmp_GM12878_TAD_filter_bed <- cbind(tmp_GM12878_TAD_filter_bed, apply(X = tmp_GM12878_TAD_filter_bed, MARGIN = 1, FUN = GetCNVofRegions, cnv_files = tmp_GM12878_CNV_40kb))
  head(tmp_GM12878_TAD_filter_bed )
  tmp_GM12878_TAD_filter_bed[, 1] <- strsplit(tmp_chr, "chr")[[1]][[2]]
  GM12878_TAD_filter_bed <- rbind(GM12878_TAD_filter_bed, tmp_GM12878_TAD_filter_bed)
  
}
# check data
dim(GM12878_TAD_filter_bed) [1] 2810    4
head(GM12878_TAD_filter_bed )
tail(GM12878_TAD_filter_bed )
write.table(x = GM12878_TAD_filter_bed, file = "GM12878_TAD_filter.bed", quote = F, row.names = F)



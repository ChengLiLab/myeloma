U266_genome_cnv_40kb <- read.table( file = "/lustre/user/liclab/wupz/DNA_seq/CNV_final/resolution_40kb/U266-merged.sorted.bam.dedup.bam_ratio.txt", 
                                    stringsAsFactors = F, header = T, sep = "\t" )
head(U266_genome_cnv_40kb)
U266_cnv_block_40kb <- FindCNVBlock( freec_ratio = U266_genome_cnv_40kb)
head(U266_cnv_block_40kb)
options(digits=14)
U266_cnv_block_40kb_bed <- cbind( chrom = paste( "chr", U266_cnv_block_40kb[, 1],sep = "" ), chromStart = as.integer(U266_cnv_block_40kb[, 2]), 
                                  chromEnd =  as.integer(U266_cnv_block_40kb[, 2] + 40000),  U266_cnv_block_40kb[, 4:6])
write.table(x = U266_cnv_block_40kb_bed, file = "/lustre/user/liclab/wupz/dosageEffect/ngs.plot/U266_cnv_block_40kb.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t")

U266_cnv_block_40kb[30, 1] == 10

U266_random_sites_40kb <- numeric()
U266_random_TAD_sites_40kb <- numeric()
for ( i in names(table(U266_cnv_block_40kb_bed[, 1]))[1:23] ) {
  tmp_length <- sum(U266_cnv_block_40kb_bed[, 1] == i)
  # generate a set of random sites used for negative control
  tmp_random_bins <- sample(x = sum(U266_all_insulation_bins[, 1] == i), size =  tmp_length )
  tmp_random_bins_sites <- U266_all_insulation_bins[U266_all_insulation_bins[, 1] == i, ][tmp_random_bins, ]
  U266_random_sites_40kb <- rbind( U266_random_sites_40kb, tmp_random_bins_sites)
  # generate a set of random TAD sites used for positive control
  tmp_random_TAD <- sample(x = sum(U266_all_tad_boundaries[, 1] == i), size =  tmp_length )
  tmp_random_TAD_sites <- U266_all_tad_boundaries[U266_all_tad_boundaries[, 1] == i, ][tmp_random_TAD, ]
  U266_random_TAD_sites_40kb <- rbind(U266_random_TAD_sites_40kb, tmp_random_TAD_sites)
}
head(U266_random_TAD_sites_40kb)
head(U266_random_sites_40kb)
write.table(x = U266_random_TAD_sites_40kb, file = "/lustre/user/liclab/wupz/dosageEffect/ngs.plot/U266_random_TAD_sites_40kb.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t")
  
write.table(x = U266_random_sites_40kb, file = "/lustre/user/liclab/wupz/dosageEffect/ngs.plot/U266_random_sites_40kb.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t" )

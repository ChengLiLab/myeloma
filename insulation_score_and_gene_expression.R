# relationship between genomic structures (insulation score) and gene expression
#1. plot gene expression along the chromosomes
ALL_exp <- read.table("/lustre/user/liclab/wupz/dosageEffect/cuffdiff_GM12878_MM/20160706_cuffdiff_hg19/genes.fpkm_tracking", header = T)
strsplit2(ALL_exp[1:10, "locus"], split = ":|-")
# find gene position
ALL_exp_gene_locus <- strsplit2(ALL_exp[, "locus"], split = ":|-")
rownames(ALL_exp_gene_locus) <- ALL_exp$gene_id
dim(ALL_exp_gene_locus)
# plot gene expression, example: chrX
plot( x = ALL_exp_gene_locus[ALL_exp_gene_locus[, 1] == "chrX", 2] , 
      y = log2(ALL_exp[ALL_exp_gene_locus[, 1] == "chrX", 14]/ALL_exp[ALL_exp_gene_locus[, 1] == "chrX", 10]) - 1,
      type = "h")
# plot CNV ration
head(ALL_exp_gene_locus)
tail(ALL_exp_gene_locus)


GM12878_is_40kb <- list()
for ( i in 1:23 ) {
  GM12878_is_40kb[[i]] <- read.table( paste('/lustre/user/liclab/lirf/Project/hic/2014/gm12878_hind3_dCTP_296mb/resolution_40k/cis/TAD_boundary/chr',
                                             i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation', sep = ''),
                                       header = T,
                                       stringsAsFactors = F)
}

RPMI8226_is_40kb <- list()
for ( i in 1:23 ) {
  RPMI8226_is_40kb[[i]] <- read.table( paste('/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release10/resolution_40k/cis/TAD_boundary/chr',
                                             i, '_chr', i, '_40k_normalmatrix.txt.is1000001.ids240001.insulation', sep = ''),
                                       header = T,
                                       stringsAsFactors = F)
}

plot(RPMI8226_is_40kb[[i]][, 8] - GM12878_is_40kb[[i]][, 8])

# example of chr2
# according to the medium of the gene, find its insulationscore in RPMI8226 and U266, GM12878


plot( x = ALL_exp_gene_locus[ALL_exp_gene_locus[, 1] == "chr2", 2] , 
      y = log2(ALL_exp[ALL_exp_gene_locus[, 1] == "chr2", 14]/ALL_exp[ALL_exp_gene_locus[, 1] == "chr2", 10]) - 1,
      type = "h", xlim = c(0, 20000000) )
lines( x = RPMI8226_is_40kb[[2]][tmp_chr2_gene_medium, 4],
       y = RPMI8226_is_40kb[[2]][, 8] - U266_is_40kb[[2]][, 8], type = "h", col = "red")

# gene expression vs insulation score in RPMI8226

tmp_chr2_gene_medium <- as.numeric(ALL_exp_gene_locus[ALL_exp_gene_locus[, 1] == "chr2", 2])/2 + as.numeric(ALL_exp_gene_locus[ALL_exp_gene_locus[, 1] == "chr2", 3])/2
tmp_insulation_RPMI8226_chr2 <- RPMI8226_is_40kb[[2]][floor(tmp_chr2_gene_medium/40000), 8]
tmp_insulation_U266_chr2 <- U266_is_40kb[[2]][floor(tmp_chr2_gene_medium/40000), 8]
tmp_insulation_GM12878_chr2 <- GM12878_is_40kb[[2]][floor(tmp_chr2_gene_medium/40000), 8]

# filtering lower expression
tmp_mask <- ALL_exp[ALL_exp_gene_locus[, 1] == "chr2", 10] == 0 | 
            is.infinite(log2(ALL_exp[ALL_exp_gene_locus[, 1] == "chr2", 14]/ALL_exp[ALL_exp_gene_locus[, 1] == "chr2", 10])) | 
            is.na(log2(ALL_exp[ALL_exp_gene_locus[, 1] == "chr2", 14]/ALL_exp[ALL_exp_gene_locus[, 1] == "chr2", 10]))

plot()

InsulationScoreChanges_RPMI8226_GM12878 <- RPMI8226_is_40kb[[2]][, 8] - GM12878_is_40kb[[2]][, 8]
InsulationScoreChanges_RPMI8226_GM12878[InsulationScoreChanges_RPMI8226_GM12878 >= 1] <- 1
InsulationScoreChanges_RPMI8226_GM12878[InsulationScoreChanges_RPMI8226_GM12878 <= -1] <- -1
plot( GM12878_is_40kb[[2]][, 2], 
      y = InsulationScoreChanges_RPMI8226_GM12878, 
      type = "h", xlab = "chromosome 2", 
      ylab = "RPMI-8226 minus GM12878",
      main = "Insulation Difference",
      col = ifelse(InsulationScoreChanges_RPMI8226_GM12878>0, "red", "blue"))

tmp_pos <- tmp_chr2_gene_medium[!tmp_mask]
tmp_is <- tmp_insulation_RPMI8226_chr2[!tmp_mask] - tmp_insulation_GM12878_chr2[!tmp_mask]
tmp_exp <- log2(ALL_exp[ALL_exp_gene_locus[, 1] == "chr2", 14]/ALL_exp[ALL_exp_gene_locus[, 1] == "chr2", 10])[!tmp_mask]
plot( x = tmp_pos, 
      y = tmp_exp, type = "h",
      xlab = "chromosome 2", 
      ylab = "log2( RPMI-8226/GM12878)",
      main = "Gene expression fold changes",
      col = ifelse(tmp_exp>0, "red", "blue")
      )

boxplot(tmp_exp[tmp_is >= 0.2], tmp_exp[tmp_is<= -0.3], 
        col = c("red", "blue"),
        main = "Relations Between Different Insulation\nand Different Expression\nT test p value: 0.0003592", 
        names = c("Different Insulation >= 0.2", "Different Insulation <= -0.3"),
        ylab = "Different Expression Fold Changes")
t.test(tmp_exp[tmp_is >= 0.2], tmp_exp[tmp_is<= -0.3]) # data:  tmp_exp[tmp_is >= 0.2] and tmp_exp[tmp_is <= -0.3] t = 3.6777, df = 115.12, p-value = 0.0003592



tmp_mask <- is.infinite(log2(ALL_exp[ALL_exp_gene_locus[, 1] == "chr2", 14]/ALL_exp[ALL_exp_gene_locus[, 1] == "chr3", 18])) | is.na(log2(ALL_exp[ALL_exp_gene_locus[, 1] == "chr3", 14]/ALL_exp[ALL_exp_gene_locus[, 1] == "chr3", 18]))
tmp_x <- tmp_insulation_RPMI8226_chr2[!tmp_mask] - tmp_insulation_U266_chr2[!tmp_mask]
tmp_y <- log2(ALL_exp[ALL_exp_gene_locus[, 1] == "chr3", 14]/ALL_exp[ALL_exp_gene_locus[, 1] == "chr3", 18])[!tmp_mask]
tmp_lm <- lm( tmp_y ~ tmp_x)
tmp_lm
summary(tmp_lm)
tmp_y2 <- ALL_exp[ALL_exp_gene_locus[, 1] == "chr3", 14][!tmp_mask]
tmp_x2 <- tmp_insulation_RPMI8226_chr2[!tmp_mask]
tmp_lm2 <- lm( tmp_y2 ~ tmp_x2)
tmp_lm2
summary(tmp_lm2)
plot( x = tmp_chr2_gene_medium[!tmp_mask] , y =  tmp_y, type = "h", xlim = c(0, 100000000))
lines( x = tmp_chr2_gene_medium[!tmp_mask] , y = tmp_x, type = "h", col = "red")

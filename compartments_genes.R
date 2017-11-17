library(diffHic)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(limma)

seg.frags <- segmentGenome(BSgenome.Hsapiens.UCSC.hg19, size=500000)
seg.frags = seg.frags[1:6087]

# read hg19 gtf file
hg19_gene_gtf <- read.table( file = "/lustre/user/liclab/publicData/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf", sep = "\t")
hg19_gene_gtf_gr <- GRanges(seqnames = hg19_gene_gtf[, 1], 
                            ranges = IRanges( start = hg19_gene_gtf[, 4], end = hg19_gene_gtf[, 5] ),
                            strand = hg19_gene_gtf[, 7] )
mcols(hg19_gene_gtf_gr ) <- DataFrame( type = hg19_gene_gtf[, 3], gene_id = strsplit2( hg19_gene_gtf[, 9], " |;")[, 2])
head(hg19_gene_gtf_gr )

# RPMI8226_HindIII
RPMI8226_HindIII_compartment <- numeric()
for(i in 1:23){
  tmp <- read.table(paste0("/lustre/user/liclab/lirf/Project/hic/mm/pca/500kb/RPMI8226_HindIII/compartment/chr",i,".txt",collapse = ""))
  RPMI8226_HindIII_compartment = rbind(RPMI8226_HindIII_compartment, tmp)
}
dim(RPMI8226_HindIII_compartment)
RPMI8226_HindIII_compartment <- as.numeric(as.matrix(RPMI8226_HindIII_compartment))
length(RPMI8226_HindIII_compartment)

# U266_HindIII
U266_HindIII_compartment <- numeric()
for(i in 1:23){
  tmp <- read.table(paste0("/lustre/user/liclab/lirf/Project/hic/mm/pca/500kb/U266_HindIII/compartment/chr",i,".txt",collapse = ""))
  U266_HindIII_compartment = rbind(U266_HindIII_compartment, tmp)
}
dim(U266_HindIII_compartment)
U266_HindIII_compartment <- as.numeric(as.matrix(U266_HindIII_compartment) )
length(U266_HindIII_compartment)

# GM12878
GM12878_compartment <- numeric()
for(i in 1:23){
  tmp <- read.table(paste0("/lustre/user/liclab/lirf/Project/hic/mm/pca/500kb/gm12878/compartment/chr",i,".txt",collapse = ""))
  GM12878_compartment = rbind(GM12878_compartment, tmp)
}
dim(GM12878_compartment)
GM12878_compartment <- as.numeric(as.matrix(GM12878_compartment))
length(GM12878_compartment)

############################################################

mcols(seg.frags) <- DataFrame(RPMI8226_HindIII = RPMI8226_HindIII_compartment,
                             U266_HindIII = U266_HindIII_compartment, GM12878 = GM12878_compartment )

seg.frags = seg.frags[which(is.na(seg.frags$RPMI8226_HindIII) == F)]
seg.frags = seg.frags[which(is.na(seg.frags$U266_HindIII) == F)]
seg.frags = seg.frags[which(is.na(seg.frags$GM12878) == F)]

seg.frags_all_A <- seg.frags[ seg.frags$RPMI8226_HindIII > 0 & seg.frags$U266_HindIII > 0 & seg.frags$GM12878 < 0]
length(seg.frags_all_A) # [1] 262

B2A_genes <- hg19_gene_gtf_gr[ findOverlaps(query = seg.frags_all_A, subject = hg19_gene_gtf_gr )@subjectHits ]

head(B2A_genes)
write.table( x = unique(B2A_genes$gene_id), file = "Compartment_B2A_in_MM.txt", quote = F, row.names = F, col.names = F)

seg.frags_all_B <- seg.frags[ seg.frags$RPMI8226_HindIII < 0 & seg.frags$U266_HindIII < 0 & seg.frags$GM12878 > 0]
length(seg.frags_all_B)

A2B_genes <- hg19_gene_gtf_gr[ findOverlaps(query = seg.frags_all_B, subject = hg19_gene_gtf_gr )@subjectHits ]
head(A2B_genes)
write.table( x = unique(A2B_genes$gene_id), file = "Compartment_A2B_in_MM.txt", quote = F, row.names = F, col.names = F) 


hg19_len <- read.table("/lustre/user/liclab/biotools/FREEC_Linux64/data/hg19.len")

# read expression of GM12878, RPMI-8226, U266
ALL_exp <- read.table("/lustre/user/liclab/wupz/dosageEffect/cuffdiff_GM12878_MM/20160706_cuffdiff_hg19/genes.fpkm_tracking", header = T)
head(ALL_exp)
rownames(ALL_exp) <- ALL_exp[,1]

plot(density(log2(as.matrix(ALL_exp[, "GM12878_FPKM"]))) )
lines( density(log2(as.matrix(ALL_exp[, "RMPI8226_FPKM"]))), col = "red")
lines( density(log2(as.matrix(ALL_exp[, "U266_FPKM"]))), col = "blue")
boxplot( log2(cbind( ALL_exp[, "GM12878_FPKM"],
                     ALL_exp[, "RMPI8226_FPKM"],
                     ALL_exp[, "U266_FPKM"] )) )


A2B_genes_exp <- ALL_exp[ unique(A2B_genes$gene_id),]
head(A2B_genes_exp)
boxplot(log2(A2B_genes_exp[, c(10,14,18)]) )
# filter the A2B genes:
# 1. GM12878 expression >= 1
# 2. log2 fold change > 2
A2B_genes_exp_filter <- A2B_genes_exp[, 10] >= 1 & log2(A2B_genes_exp[, 10]/A2B_genes_exp[, 14]) >= 1 & log2(A2B_genes_exp[, 10]/A2B_genes_exp[, 18]) >= 1
sum(A2B_genes_exp_filter)

write.table( rownames( A2B_genes_exp[ A2B_genes_exp_filter, ] ), 
             file = "Compartment_A2B_in_MM_diff_exp.txt", quote = F, row.names = F, col.names = F )

B2A_genes_exp <- ALL_exp[ intersect(unique(B2A_genes$gene_id), ALL_exp[,1]), ]
boxplot(log2(B2A_genes_exp[, c(10,14,18)]) )

head(B2A_genes_exp)
dim(B2A_genes_exp)
# filter the A2B genes:
# 1. GM12878 expression >= 1
# 2. log2 fold change > 2
B2A_genes_exp_filter <- B2A_genes_exp[, 14] >= 1 & B2A_genes_exp[, 18] >= 1 & log2(B2A_genes_exp[, 14]/B2A_genes_exp[, 10]) >= 1 & log2(B2A_genes_exp[, 18]/B2A_genes_exp[, 10]) >= 1
sum(B2A_genes_exp_filter)
B2A_genes_exp[ B2A_genes_exp_filter, ]
write.table( rownames( B2A_genes_exp[ B2A_genes_exp_filter, ] ), 
             file = "Compartment_B2A_in_MM_diff_exp.txt", quote = F, row.names = F, col.names = F )

ALL_exp[ALL_exp[, 1] == "MAP4K3", ]
# find the location of genes
protein_coding_gene <- read.csv("/lustre/user/liclab/wupz/database/human_genes/protein-coding_gene.csv")
a2b_enriched_genes <- c( "SPRY2", "CBLB", "PIK3R5", "IL7R", 
                         "IL1R2", "IL1R1", "IL7R", 
                         "IL18R1", "IL1R2", "IL1R1", "IL7R" )
a2b_enriched_genes_matrix <- numeric()
for (i in b2a_enriched_genes) {
  a2b_enriched_genes_matrix <- rbind(a2b_enriched_genes_matrix, protein_coding_gene[protein_coding_gene[, 2] == i, ])
}
write.csv(x = a2b_enriched_genes_matrix, file = "a2b_enriched_genes_matrix.txt", 
            quote = F, row.names = F )


b2a_enriched_genes <- c( "MAP4K3", 
                         "CACNG6", 
                         "FGF13", 
                         "DUSP6", 
                         "NCAM1", 
                         "CADM1", 
                         "NLGN1", 
                         "IRS2", 
                         "PDE3B", 
                         "PYGB" )
b2a_enriched_genes_matrix <- numeric()
for (i in b2a_enriched_genes) {
  b2a_enriched_genes_matrix <- rbind(b2a_enriched_genes_matrix, protein_coding_gene[protein_coding_gene[, 2] == i, ])
}
write.csv(x = b2a_enriched_genes_matrix, file = "b2a_enriched_genes_matrix.txt", 
          quote = F, row.names = F )


MM_genes <- read.table("../MM_genes/mm_related_keygene_in_RNA_seq.txt", header = T)
MM_genes_a2b_overlap <- intersect(x = as.vector(MM_genes[, 6]), y = unique(A2B_genes$gene_id) )
MM_genes_b2a_overlap <- intersect(x = as.vector(MM_genes[, 6]), y = unique(B2A_genes$gene_id) )
ALL_exp[as.vector(MM_genes[, 6]), c(7, 10,14,18)]
MM_genes_gr <- GRanges(seqnames = MM_genes[, "chr"], 
                       ranges = IRanges( start = MM_genes[, "start"], 
                                         end = MM_genes[, "end"]),
                       strand = MM_genes[, "strand"] )
mcols(MM_genes_gr ) <- DataFrame( gene = MM_genes[, "gene"], gene_id = MM_genes[, "gene_id"])
MM_genes_gr

mcols(MM_genes_gr)$RPMI8226 <- 0
mcols(MM_genes_gr)$U266 <- 0
mcols(MM_genes_gr)$GM12878 <- 0
mcols(MM_genes_gr)[findOverlaps(query = promoters(MM_genes_gr), subject = seg.frags)@queryHits, 3:5] <- mcols(seg.frags[findOverlaps(query = promoters(MM_genes_gr), subject = seg.frags)@subjectHits])


sum(MM_genes_gr$RPMI8226 > 0 & MM_genes_gr$U266 > 0 & MM_genes_gr$GM12878 > 0) # 27
sum(MM_genes_gr$RPMI8226 > 0 & MM_genes_gr$U266 < 0 & MM_genes_gr$GM12878 > 0) # 1
sum(MM_genes_gr$RPMI8226 < 0 & MM_genes_gr$U266 > 0 & MM_genes_gr$GM12878 > 0) # 5
sum(MM_genes_gr$RPMI8226 < 0 & MM_genes_gr$U266 < 0 & MM_genes_gr$GM12878 > 0) # 0
sum(MM_genes_gr$RPMI8226 < 0 & MM_genes_gr$U266 < 0 & MM_genes_gr$GM12878 < 0) # 8
sum(MM_genes_gr$RPMI8226 < 0 & MM_genes_gr$U266 > 0 & MM_genes_gr$GM12878 < 0) # 4
sum(MM_genes_gr$RPMI8226 > 0 & MM_genes_gr$U266 < 0 & MM_genes_gr$GM12878 < 0) # 0
sum(MM_genes_gr$RPMI8226 > 0 & MM_genes_gr$U266 > 0 & MM_genes_gr$GM12878 < 0) # 0


enrichment_results_of_all_genes_with_compartment_switching <- read.table(file = "../comparments_genes/enrichment_results_of_all_genes_with_compartment_switching.txt", sep = "\t", header = T)
head(enrichment_results_of_all_genes_with_compartment_switching)
dim(enrichment_results_of_all_genes_with_compartment_switching)

png(filename = "../comparments_genes/enrichment_results_of_all_genes_with_compartment_switching.png")
par(cex = 2)
barplot( -log10(enrichment_results_of_all_genes_with_compartment_switching[, "PValue"]), 
        horiz = T, col = "darkorange", xlab = "-log10 (p-value)")
dev.off()

# 
enrichment_results_of_all_genes_with_compartment_switching_and_diff_exp <- read.table(file = "../comparments_genes/enrichment_results_of_all_genes_with_compartment_switching_and_diff_exp.txt", sep = "\t", header = T)
head(enrichment_results_of_all_genes_with_compartment_switching_and_diff_exp)
dim(enrichment_results_of_all_genes_with_compartment_switching_and_diff_exp)

png(filename = "../comparments_genes/enrichment_results_of_all_genes_with_compartment_switching_and_diff_exp.png")
par(cex = 2)
barplot( -log10(enrichment_results_of_all_genes_with_compartment_switching[1:4, "PValue"]), 
         horiz = T, col = "darkorange", xlab = "-log10 (p-value)")
dev.off()

enrichment_results_A2B_enrichr_KEGG <- read.table(file = "../comparments_genes/A2B_KEGG_2016_table.txt", sep = "\t", header = T)
enrichment_results_B2A_enrichr_KEGG <- read.table(file = "../comparments_genes/B2A_KEGG_2016_table.txt", sep = "\t", header = T)
sum(enrichment_results_A2B_enrichr_KEGG[, 3]  <= 0.05) # 21
sum( enrichment_results_B2A_enrichr_KEGG[, 3] <= 0.05 ) # 7
enrichment_results_enrichr_KEGG <- rbind( enrichment_results_A2B_enrichr_KEGG[1:21, c(1,3)], enrichment_results_B2A_enrichr_KEGG[1:7, c(1,3)])
library(limma)
enrichment_results_enrichr_KEGG[,1] <- as.vector( strsplit2(enrichment_results_enrichr_KEGG[,1], "_")[,1])

barplot(height = enrichment_results_enrichr_KEGG[, 2], names.arg = T)
head(enrichment_results_B2A_enrichr_KEGG)
head(enrichment_results_A2B_enrichr_KEGG )
dim(enrichment_results_of_all_genes_with_compartment_switching)

png(filename = "../comparments_genes/enrichment_results_of_all_genes_with_compartment_switching.png", width = 2048, height = 1024)
par(cex = 2, mar= c(5, 30, 3, 3))
barplot( -log10(enrichment_results_enrichr_KEGG[, 2]), 
         names.arg = enrichment_results_enrichr_KEGG[,1],las=1,
         horiz = T,  xlab = "-log10 (p-value)",
         col = c( rep("deepskyblue1", 21), rep("orangered1", 7)) )
dev.off()



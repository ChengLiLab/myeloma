# 20160629
# wupz
# read gtf file and find genes in specific regions


# 1. read gtf
hg19_gene_gtf <- read.table( file = "/lustre/user/liclab/publicData/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf", sep = "\t")

findGene <- function ( query, gene_gtf) {
  chri <- query[1]
  chri_gene_gtf <- gene_gtf[gene_gtf[, 1] == chri, ]
  find_id <- (chri_gene_gtf[, 4] >= query[2]) & (chri_gene_gtf[, 5] <= query[3])
  chri_gene_gtf_find <- chri_gene_gtf[find_id, ]
  return(chri_gene_gtf_find)
}

head(hg19_gene_gtf)
# test
test_region <- c("chr11", 47380528, 47380528)
test_find <- findGene( query = test_region, gene_gtf = hg19_gene_gtf)
test_find
dim(test_find)

RPMI_CTX_validate <- RPMI8226_CTX_filtered_sort_dedup_2[ctx_score[, 5] >= 20, ]
class(RPMI_CTX_validate)
RPMI_CTX_validate_region1 <- cbind( RPMI_CTX_validate[, 1], 
                                    as.numeric(RPMI_CTX_validate[, 2]) - 50000,  
                                    as.numeric(RPMI_CTX_validate[, 2]) + 50000 )
RPMI_CTX_validate_region1
RPMI_CTX_validate_region1_genes <- list()
for ( i in 1:dim(RPMI_CTX_validate_region1)[1] ) {
  RPMI_CTX_validate_region1_genes[[i]] <- list()
  RPMI_CTX_validate_region1_genes[[i]][[1]] <- RPMI_CTX_validate[i, ]
  RPMI_CTX_validate_region1_genes[[i]][[2]] <- findGene(RPMI_CTX_validate_region1[i, ], gene_gtf = hg19_gene_gtf)
}

RPMI_CTX_validate_region1_genes <- findGene(RPMI_CTX_validate_region1[1, ], gene_gtf = hg19_gene_gtf)
dim(RPMI_CTX_validate_region1_genes)

RPMI_CTX_validate_region2 <- cbind( RPMI_CTX_validate[, 5], 
                                    floor(as.numeric(RPMI_CTX_validate[, 6])/40000) * 40000,  
                                    ceiling(as.numeric(RPMI_CTX_validate[, 6])/40000) * 40000 )
RPMI_CTX_validate_region2
RPMI_CTX_validate_region2_genes <- findGene(RPMI_CTX_validate_region2[1, ], gene_gtf = hg19_gene_gtf)
dim(RPMI_CTX_validate_region2_genes)

RPMI_CTX_validate_region2_genes <- list()
for ( i in 1:dim(RPMI_CTX_validate_region2)[1] ) {
  RPMI_CTX_validate_region2_genes[[i]] <- list()
  RPMI_CTX_validate_region2_genes[[i]][[1]] <- RPMI_CTX_validate[i, ]
  RPMI_CTX_validate_region2_genes[[i]][[2]] <- findGene(RPMI_CTX_validate_region2[i, ], gene_gtf = hg19_gene_gtf)
}
RPMI_CTX_validate_region2_genes[[3]]

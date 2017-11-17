RPMI8226_CTX <- read.table(file = "/lustre/user/liclab/wupz/P3_HiC/translocation/RPMI8226/RPMI8226_CTX.txt", 
                           stringsAsFactors = F )
dim(RPMI8226_CTX)
apply(RPMI8226_CTX, 2, class)

head(RPMI8226_CTX)

plot(density(RPMI8226_CTX[, 8]/RPMI8226_CTX[, 11]) )

tmp_filter <- RPMI8226_CTX[, 8]/RPMI8226_CTX[, 11] >= 0.1 & RPMI8226_CTX[, 4]/RPMI8226_CTX[, 10] >= 0.1
sum(tmp_filter)

RPMI8226_CTX_filtered <- RPMI8226_CTX[tmp_filter, ]
table(RPMI8226_CTX_filtered[,1])
tmp_filter2 <- RPMI8226_CTX_filtered[,1] != "chrX" & RPMI8226_CTX_filtered[,1] != "chrY" & RPMI8226_CTX_filtered[,5] != "chrX" & RPMI8226_CTX_filtered[,5] != "chrY"
sum(tmp_filter2)
RPMI8226_CTX_filtered_2 <- RPMI8226_CTX_filtered[tmp_filter2, ]

getScore <- function(query, hic_matrix) {
  chr1 <- as.integer(strsplit(query[1], split = "chr")[[1]][2])
  bin1 <- as.integer( query[2] / 40000)
  chr2 <- as.integer(strsplit(query[5], split = "chr")[[1]][2])
  bin2 <- as.integer( query[6] / 40000)
  if ( chr2 > chr1) {
    tmp  <- chr2
    chr2 <- chr1
    chr1 <- tmp
    tmp <- bin2
    bin2 <- bin1
    bin1 <- bin2
  }
  score <- hic_matrix[[]][[]]
}
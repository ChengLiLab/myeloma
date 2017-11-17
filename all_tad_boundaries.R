# 20160424
# Read TAD boundaries from 23 chromosomes and combine to one file
# wupz

all_tad_boundaries <- numeric()
for ( i in 1:23) {
  tmp_tad_boundaries <- read.table( paste("/lustre/user/liclab/lirf/Project/hic/data.2015.6.24/release12.2/resolution_40k/cis/TAD_boundary/chr",
                                          i, "_chr", i, "_40k_normalmatrix.txt.is1000001.ids240001.insulation.boundaries.bed", sep = ""),
                                    skip = 1
    )
  head(tmp_tad_boundaries)
  all_tad_boundaries <- rbind(all_tad_boundaries, tmp_tad_boundaries)
}

write.table( x = all_tad_boundaries, file = "/lustre/user/liclab/wupz/dosageEffect/ngs.plot/U266_HindIII_40kb_TAD_boundary.bed", 
             sep = "\t", row.names = F, quote = F, col.names = F)

dim(all_tad_boundaries)
tmp_tad_414 <- sample(4112, 434)
write.table( x = all_tad_boundaries[tmp_tad_414, ], file = "/lustre/user/liclab/wupz/dosageEffect/ngs.plot/U266_HindIII_40kb_TAD_boundary_434_sites.bed", 
             sep = "\t", row.names = F, quote = F, col.names = F)

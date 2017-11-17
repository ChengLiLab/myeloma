## 3.6 the relationship between TADs and A/B 
## WUPZ
## 20170511

## 1. compartment A/B
### 1.1 Read data
U266_compartment_MobI <- vector()
for (i in 1:23) {
  tmp_U266_compartment_MobI <- read.table(paste0("/lustre/user/liclab/lirf/Project/hic/mm/pca/500kb/U266_MboI/compartment/chr",i,".txt"),
                                          stringsAsFactors =  F )
  tmp_length <- dim(tmp_U266_compartment_MobI)[1]
  tmp_start <- 0:(tmp_length-1) * 500000
  tmp_end <- 1:tmp_length * 500000
  tmp_U266_compartment_MobI_bed <- cbind( chr = i,
                                      start = tmp_start,
                                      end = tmp_end,
                                      pc1 = tmp_U266_compartment_MobI
                                      )
  U266_compartment_MobI <- rbind(U266_compartment_MobI, tmp_U266_compartment_MobI_bed )
}

U266_compartment_MobI <- cbind(U266_compartment_MobI, compartment = NA )
U266_compartment_MobI[which(U266_compartment_MobI[, 4] > 0), "compartment"] <- "A"
U266_compartment_MobI[which(U266_compartment_MobI[, 4] < 0), "compartment"] <- "B"
### check data
dim(U266_compartment_MobI)
head(U266_compartment_MobI)
### 1.2 Count compartment A/B
U266_compartment_MobI_bed <- FindCompartment(U266_compartment_MobI)
png(filename = "Distribution of length of compartment AB in U266.png", width = 1024, height = 1024)
tmp_compartment_length <- U266_compartment_MobI_bed$end-U266_compartment_MobI_bed$start
par(cex = 3, lwd = 3, mgp = c(2.5,1,0))
plot(density(tmp_compartment_length), 
     main = "Distribution of compartment length, U266", 
     xlab = paste("N = ", length(tmp_compartment_length), ", Median length: ", median(tmp_compartment_length), "bp"))
dev.off()
## 2. Define A/B boudary
### if there is a A-B switch, define the 40 kb border of the two compartment as their boudary
tmp_boundary_id <- 1
U266_MobI_compartment_boundaries_bed <- vector()
for ( i in 2:dim(U266_compartment_MobI_bed)[1] ) {
  print(tmp_boundary_id)
  if( U266_compartment_MobI_bed[i, 5] != U266_compartment_MobI_bed[i-1, 5] & (U266_compartment_MobI_bed[i, 1] == U266_compartment_MobI_bed[i-1, 1]) ) {
    tmp_U266_MobI_compartment_boundaries_bed <- cbind(chr = U266_compartment_MobI_bed[i, 1],
                                                      start = U266_compartment_MobI_bed[i, 2] - 20000,
                                                      end = U266_compartment_MobI_bed[i, 2] + 20000,
                                                      name = paste0("Boundary_", tmp_boundary_id)
                                                     )
    U266_MobI_compartment_boundaries_bed <- rbind(U266_MobI_compartment_boundaries_bed, tmp_U266_MobI_compartment_boundaries_bed)
    tmp_boundary_id <- tmp_boundary_id + 1
  }
}
write.table(x = U266_MobI_compartment_boundaries_bed, 
            file = "U266_MobI_compartment_boundaries.bed", 
            quote = F, 
            row.names = F,
            sep = "\t", col.names = F)

## 3. Test if TAD boundaries are overlaped with compartment A/B
### read distance data
U266_TAD_boundary_AB_boundary_distance_bed <- read.table("U266_TAD_boundary_AB_boundary_distance.bed")
head(U266_TAD_boundary_AB_boundary_distance_bed)
tmp_distance <- floor(U266_TAD_boundary_AB_boundary_distance_bed[, 10]/40000)
table(tmp_distance)
png(filename = "Distance of AB compartment boundaries and TAD boundaries in U266.png", width = 1024, height = 1024)
par(cex = 3, lwd = 3, mgp = c(2.5,1,0))
plot(density(tmp_distance),
     main = "Distance of AB compartment boundaries\n and TAD boundaries, U266", 
     xlab = paste("N = ", length(tmp_distance), ", Median length: ", median(tmp_distance), "bin"))
dev.off()


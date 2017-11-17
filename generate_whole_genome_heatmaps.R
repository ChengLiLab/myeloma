# generate whole genome raw matrix with 5M bin resolution
# 20151220
# wupz
# the raw data is 200kb resolution
# first change the resolution to 5M,
# then combine all cis and trans interaction matrix to a whole matrix
# third draw the heatmap
library(ComplexHeatmap)
library(circlize)
source('/lustre/user/liclab/wupz/dosageEffect/scripts/ChangeResolution.R')
source('/lustre/user/liclab/wupz/dosageEffect/scripts/GenerateWholeGenomeMatrix.R')
# 1.1 change resolution to 5M
raw_HiC_matrix_5M <- list()
for (i in 1:length(raw_HiC_matrix) ) {
  raw_HiC_matrix_5M[[i]] <- list()
  for (j in 1:length(raw_HiC_matrix[[i]])) {
    raw_HiC_matrix_5M[[i]][[j]] <- list()
    for ( k in j:length(raw_HiC_matrix[[i]][[j]])) {
      if( k == j ) {
        raw_HiC_matrix_5M[[i]][[j]][[k]] <- ChangeResolution(raw_HiC_matrix[[i]][[j]][[k]], 
                                                             resolution_change_folds = 25 , 
                                                             cis = T)
      }
      else 
        raw_HiC_matrix_5M[[i]][[j]][[k]] <- ChangeResolution(raw_HiC_matrix[[i]][[j]][[k]], 
                                                             resolution_change_folds = 25 , 
                                                             cis = F)
      
    }
  }
}
save(raw_HiC_matrix_5M, file = '/lustre/user/liclab/wupz/dosageEffect/preprocessingData/raw_HiC_matrix_5M.Rdata')

# 1.2 change resolution to 2M
raw_HiC_matrix_2M <- list()
for (i in 1:length(raw_HiC_matrix) ) {
  raw_HiC_matrix_2M[[i]] <- list()
  for (j in 1:length(raw_HiC_matrix[[i]])) {
    raw_HiC_matrix_2M[[i]][[j]] <- list()
    for ( k in j:length(raw_HiC_matrix[[i]][[j]])) {
      if( k == j ) {
        raw_HiC_matrix_2M[[i]][[j]][[k]] <- ChangeResolution(raw_HiC_matrix[[i]][[j]][[k]], 
                                                             resolution_change_folds = 10 , 
                                                             cis = T)
      }
      else 
        raw_HiC_matrix_2M[[i]][[j]][[k]] <- ChangeResolution(raw_HiC_matrix[[i]][[j]][[k]], 
                                                             resolution_change_folds = 10 , 
                                                             cis = F)
      
    }
  }
}
save(raw_HiC_matrix_2M, file = '/lustre/user/liclab/wupz/dosageEffect/preprocessingData/raw_HiC_matrix_2M.Rdata')

# 1.2 change resolution to 2M
raw_HiC_matrix_1M <- list()
for (i in 1:length(raw_HiC_matrix) ) {
  raw_HiC_matrix_1M[[i]] <- list()
  for (j in 1:length(raw_HiC_matrix[[i]])) {
    raw_HiC_matrix_1M[[i]][[j]] <- list()
    for ( k in j:length(raw_HiC_matrix[[i]][[j]])) {
      if( k == j ) {
        raw_HiC_matrix_1M[[i]][[j]][[k]] <- ChangeResolution(raw_HiC_matrix[[i]][[j]][[k]], 
                                                             resolution_change_folds = 5 , 
                                                             cis = T)
      }
      else 
        raw_HiC_matrix_1M[[i]][[j]][[k]] <- ChangeResolution(raw_HiC_matrix[[i]][[j]][[k]], 
                                                             resolution_change_folds = 5 , 
                                                             cis = F)
      
    }
  }
}
save(raw_HiC_matrix_1M, file = '/lustre/user/liclab/wupz/dosageEffect/preprocessingData/raw_HiC_matrix_1M.Rdata')
RPMI8226_HindIII_1M <- GenerateWholeGenomeMatrix (  raw_HiC_matrix_1M[[1]] )
RPMI8226_MboI_1M <- GenerateWholeGenomeMatrix (  raw_HiC_matrix_1M[[2]] )
U266_HindIII_1M <- GenerateWholeGenomeMatrix (  raw_HiC_matrix_1M[[4]] )
U266_MboI_1M <- GenerateWholeGenomeMatrix (  raw_HiC_matrix_1M[[3]] )
write.table(x = RPMI8226_HindIII_1M[[2]] , file = "/lustre/user/liclab/wupz/dosageEffect/hicNormaliztion/1000kb/RPMI8226_HindIII_1000kb.txt", 
            row.names = F, col.names = F, sep="\t", quote=F )
write.table(x = RPMI8226_MboI_1M[[2]], file = "/lustre/user/liclab/wupz/dosageEffect/hicNormaliztion/1000kb/RPMI8226_MboI_1000kb.txt", 
            row.names = F, col.names = F, sep="\t", quote=F )
write.table(x = U266_MboI_1M[[2]], file = "/lustre/user/liclab/wupz/dosageEffect/hicNormaliztion/1000kb/U266_MboI_1000kb.txt", 
            row.names = F, col.names = F, sep="\t", quote=F )
write.table(x = U266_HindIII_1M[[2]], file = "/lustre/user/liclab/wupz/dosageEffect/hicNormaliztion/1000kb/U266_HindIII_1000kb.txt", 
            row.names = F, col.names = F, sep="\t", quote=F )

# 2. generate whole genome interaction matrix in 5M resolution
RPMI8226_HindIII_5M <- GenerateWholeGenomeMatrix (  raw_HiC_matrix_5M[[1]] )
RPMI8226_MboI_5M <- GenerateWholeGenomeMatrix (  raw_HiC_matrix_5M[[2]] )
U266_HindIII_5M <- GenerateWholeGenomeMatrix (  raw_HiC_matrix_5M[[4]] )
U266_MboI_5M <- GenerateWholeGenomeMatrix (  raw_HiC_matrix_5M[[3]] )

write.table(RPMI8226_HindIII_5M[[2]], file = '../heatmaps_of_whole_genome/RPMI8226_HindIII_5M.txt', col.names = F, row.names = F)
write.table(RPMI8226_MboI_5M[[2]], file = '../heatmaps_of_whole_genome/RPMI8226_MboI_5M.txt', col.names = F, row.names = F)
write.table(U266_HindIII_5M[[2]], file = '../heatmaps_of_whole_genome/U266_HindIII_5M.txt', col.names = F, row.names = F)
write.table(U266_MboI_5M[[2]], file = '../heatmaps_of_whole_genome/U266_MboI_5M.txt', col.names = F, row.names = F)

# 2.2 found the centromeres or other regions that have no reads
mask_region_vector_5M <- vector(mode = 'logical', length = 0)
mask_region_list_5M <- list()
for ( i in 1:24 ) {
  mask_region_list_5M[[i]] <- apply(raw_HiC_matrix_5M[[1]][[i]][[i]], 1, sum ) == 0 | 
    apply(raw_HiC_matrix_5M[[2]][[i]][[i]], 1, sum ) == 0 | 
    apply(raw_HiC_matrix_5M[[3]][[i]][[i]], 1, sum ) == 0 | 
    apply(raw_HiC_matrix_5M[[4]][[i]][[i]], 1, sum ) == 0 
  mask_region_vector_5M <- c(mask_region_vector_5M, mask_region_list_5M[[i]])
}
# 2.3 set the value in region that have no read to NA
RPMI8226_HindIII_5M[[2]][mask_region_vector_5M, ] <- NA
RPMI8226_HindIII_5M[[2]][, mask_region_vector_5M] <- NA
RPMI8226_MboI_5M[[2]][mask_region_vector_5M, ] <- NA
RPMI8226_MboI_5M[[2]][, mask_region_vector_5M] <- NA
U266_HindIII_5M[[2]][mask_region_vector_5M, ] <- NA
U266_HindIII_5M[[2]][, mask_region_vector_5M] <- NA
U266_MboI_5M[[2]][mask_region_vector_5M, ] <- NA
U266_MboI_5M[[2]][, mask_region_vector_5M] <- NA
# 3 draw heatmap
figure_dir <- '/lustre/user/liclab/wupz/dosageEffect/heatmaps_of_whole_genome/'
# 3.1 RPMI-8226, HindIII
png(filename = paste(figure_dir, 'RPMI8226_HindIII_5M.png', sep = '/'), width = 2048, height = 2048)
Heatmap(RPMI8226_HindIII_5M[[2]], 
        name = 'Raw\nCounts\n', 
        column_title = 'RPMI-8226, HindIII\n',
        column_title_gp = gpar(fontsize = 100), 
        col = colorRamp2( quantile(RPMI8226_HindIII_5M[[2]], c(0.05, 0.98), na.rm = T), c( "white", "red") ),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 75), 
                                    labels_gp = gpar(fontsize = 75), 
                                    color_bar = "continuous",
                                    grid_width = unit(20, "mm") )
)
dev.off()
# 3.2 RPMI-8226, MboI
png(filename = paste(figure_dir, 'RPMI8226_MboI_5M.png', sep = '/'), width = 2048, height = 2048)
Heatmap(RPMI8226_MboI_5M[[2]], 
        name = 'Raw\nCounts\n', 
        column_title = 'RPMI-8226, MboI\n',
        column_title_gp = gpar(fontsize = 100), 
        col = colorRamp2( quantile(RPMI8226_MboI_5M[[2]], c(0.05, 0.98), na.rm = T), c( "white", "red") ),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 75), 
                                    labels_gp = gpar(fontsize = 75), 
                                    color_bar = "continuous",
                                    grid_width = unit(20, "mm") )
)
dev.off()

# 3.3 RPMI-8226, MboI + HindIII
png(filename = paste(figure_dir, 'RPMI8226_5M.png', sep = '/'), width = 2048, height = 2048)
Heatmap(RPMI8226_HindIII_5M[[2]] + RPMI8226_MboI_5M[[2]], 
        name = 'Raw\nCounts\n', 
        column_title = 'RPMI 8226',
        column_title_gp = gpar(fontsize = 100), 
        col = colorRamp2( quantile(RPMI8226_HindIII_5M[[2]] + RPMI8226_MboI_5M[[2]], c(0.05, 0.98), na.rm = T), c( "white", "red") ),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 75),
                                    labels_gp = gpar(fontsize = 75),
                                    color_bar = "continuous",
                                    grid_width = unit(20, "mm") )
)
dev.off()

# 3.4 U266, HindIII
png(filename = paste(figure_dir, 'U266_HindIII_5M.png', sep = '/'), width = 2048, height = 2048)
Heatmap(U266_HindIII_5M[[2]], 
        name = 'Raw\nCounts\n', 
        column_title = 'U266, HindIII\n',
        column_title_gp = gpar(fontsize = 100), 
        col = colorRamp2( quantile(U266_HindIII_5M[[2]], c(0.05, 0.98), na.rm = T), c( "white", "red") ),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 75), 
                                    labels_gp = gpar(fontsize = 75), 
                                    color_bar = "continuous",
                                    grid_width = unit(20, "mm") )
)
dev.off()
# 3.5 U266, MboI
png(filename = paste(figure_dir, 'U266_MboI_5M.png', sep = '/'), width = 2048, height = 2048)
Heatmap(U266_MboI_5M[[2]], 
        name = 'Raw\nCounts\n', 
        column_title = 'U266, MboI\n',
        column_title_gp = gpar(fontsize = 100), 
        col = colorRamp2( quantile(U266_MboI_5M[[2]], c(0.05, 0.98), na.rm = T), c( "white", "red") ),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 75),
                                    labels_gp = gpar(fontsize = 75),
                                    color_bar = "continuous",
                                    grid_width = unit(20, "mm") )
)
dev.off()

# 3.6 U266, MboI + HindIII
png(filename = paste(figure_dir, 'U266_5M.png', sep = '/'), width = 2048, height = 2048)
Heatmap(U266_HindIII_5M[[2]] + U266_MboI_5M[[2]], 
        name = 'Raw\nCounts\n', 
        column_title = 'U266\n',
        column_title_gp = gpar(fontsize = 100), 
        col = colorRamp2( quantile(U266_HindIII_5M[[2]] + U266_MboI_5M[[2]], c(0.05, 0.98), na.rm = T), c( "white", "red") ),
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 75),
                                    labels_gp = gpar(fontsize = 75),
                                    color_bar = "continuous",
                                    grid_width = unit(20, "mm") )
)
dev.off()

# proportion of cis and trans interction
RPMI8226_cis_total_reads <- 0
RPMI8226_trans_total_reads <- 0
U266_cis_total_reads <- 0
U266_trans_total_reads <- 0
for( i in 1:24 ) {
  RPMI8226_cis_total_reads <- RPMI8226_cis_total_reads + UpperTriSum(raw_HiC_matrix_5M[[1]][[i]][[i]]) + UpperTriSum(raw_HiC_matrix_5M[[2]][[i]][[i]])
  U266_cis_total_reads <- U266_cis_total_reads + UpperTriSum(raw_HiC_matrix_5M[[3]][[i]][[i]]) + UpperTriSum(raw_HiC_matrix_5M[[4]][[i]][[i]])
}
RPMI8226_trans_total_reads <- UpperTriSum( (RPMI8226_HindIII_5M[[2]] + RPMI8226_MboI_5M[[2]])[!mask_region_vector_5M, !mask_region_vector_5M ] ) - RPMI8226_cis_total_reads
U266_trans_total_reads <- UpperTriSum( (U266_HindIII_5M[[2]] + U266_MboI_5M[[2]])[!mask_region_vector_5M, !mask_region_vector_5M ] ) - U266_cis_total_reads

png( filename = 'RPMI8226 Proportion of cis and trans.png',width = 1024, height = 1024)
par(cex = 3)
RPMI8226_cis_trans_prop <- c(RPMI8226_cis_total_reads,RPMI8226_trans_total_reads )/(RPMI8226_cis_total_reads + RPMI8226_trans_total_reads )
RPMI8226_cis_trans_prop <- paste(round(RPMI8226_cis_trans_prop, 3) *100, '%', sep = '')
pie(c(RPMI8226_cis_total_reads,RPMI8226_trans_total_reads ), 
    main = 'RPMI8226\nProportion of Cis and Trans Interactions',
    labels = paste(c('Cis', 'Trans'), RPMI8226_cis_trans_prop), 
    col = c('orangered', 'royalblue1') )
dev.off()

png( filename = 'U266 Proportion of cis and trans.png',width = 1024, height = 1024)
par(cex = 3)
U266_cis_trans_prop <- c(U266_cis_total_reads,U266_trans_total_reads )/(U266_cis_total_reads + U266_trans_total_reads )
U266_cis_trans_prop <- paste(round(U266_cis_trans_prop, 3) *100, '%', sep = '')
pie(c(U266_cis_total_reads,U266_trans_total_reads ), 
    main = 'U266\nProportion of Cis and Trans Interactions',
    labels = paste(c('Cis', 'Trans'), U266_cis_trans_prop), 
    col = c('orangered', 'royalblue1') )
dev.off()
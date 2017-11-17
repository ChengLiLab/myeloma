# Draw the heatmap of upper triangle

tri_heatmap <- function( raw_matrix,... ) {
  
  row_num <- dim(raw_matrix)[1]
  col_num <- dim(raw_matrix)[2]
  raw_matrix[is.na(raw_matrix)] <- 0
  plot(x=0:col_num, y=0:(col_num)*1/2,  ylab = "",xlab = "",
       xaxt = "n", yaxt = "n", bty = "n", frame.plot = F, type = "n", ...)
  
  for ( i in 0:(row_num*1/2) ) {
    x1 <- seq(0 + i/2, row_num-1-i/2, 1)
    x2 <- x1 + 1
    y1 <- i
    y2 <- y1 + 1 
    #生成坐标vector
    #生成颜色vector
    
    max_red  <- quantile(raw_matrix[raw_matrix!=0],probs = 0.95, na.rm = T)
    raw_vector <- diag( raw_matrix[ (i+1):row_num, 1:(col_num - i) ] )
    color_vector <- raw_vector / max_red
    color_vector[color_vector > 1 ] <- 1
    color_vector[color_vector < 0 ] <- 0
    red_vector <- rep(1,length(raw_vector))
    green_vector <- 1 - color_vector
    blue_vector <- 1 - color_vector
    rect(x1,y1,x2,y2,col = rgb(red_vector,green_vector,blue_vector), border = NA)
  }
  
}

#tri_heatmap( raw_matrix = raw_HiC_matrix[[1]][[1]][[1]][1:10,1:10])

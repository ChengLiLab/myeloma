# heatmap based on rect
#  用基本函数画热图
rect_heatmap <- function(raw_matrix, col = NULL, ... ) {
  row_num <- dim(raw_matrix)[1]
  col_num <- dim(raw_matrix)[2]
  raw_matrix[is.na(raw_matrix)] <- 0
  #生成坐标vector
  x1 <- rep(c(0:(col_num-1)),row_num)
  x2 <- x1 + 1
  y1 <- rep(c(0:(-col_num+1)),each=row_num)
  y2 <- y1 -1 
  
  plot(x=0:col_num, y=-row_num:0, type="n", frame.plot = F, asp = 1, 
       xaxt = "n", yaxt = "n", ... )
  #生成颜色vector
  if ( is.null(col) ) {
    region <- quantile(raw_matrix[!(raw_matrix == 0)],probs = c(0,0.95,1), na.rm = T)
    raw_vector <- as.vector( t(raw_matrix) )
    color_vector <- raw_vector / floor(region[2])
    color_vector[color_vector > 1 ] <- 1
    red_vector <- rep(1,length(raw_vector))
    green_vector <- 1 - color_vector
    blue_vector <- 1 - color_vector
  }
  rect(x1,y1,x2,y2,col = rgb(red_vector,green_vector,blue_vector), border = NA)
}

add_legend <- function ( raw_matrix, ...) {
  raw_matrix[is.na(raw_matrix)] <- 0
  tmp_raw_matrix_vector <- as.vector(as.matrix(raw_matrix) )
  # set 0.95 as the max red color scale
  tmp_raw_matrix_vector_0.95 <- quantile(tmp_raw_matrix_vector[!(tmp_raw_matrix_vector == 0)], 0.95 )
  #tmp_raw_matrix_vector_NonZero <- tmp_raw_matrix_vector[tmp_raw_matrix_vector >0 & tmp_raw_matrix_vector <= tmp_raw_matrix_vector_0.95]
  # calculate 5 red colors
  legend_value <- ceiling( tmp_raw_matrix_vector_0.95 * c(1:5)/5 )
  # par(mar=c(10,0,10,2))
  plot(x= 1:10,y= seq(0,20, length.out = 10),type="n",frame.plot = F,ylab = "",xlab = "",yaxt = "n",xaxt = "n", ...)
  legend_x1 <- rep(2.5, 5) 
  legend_x2 <- rep(5, 5)
  legend_y1 <- seq(8, 12,length.out =5)
  legend_y2 <- seq(8,12,length.out =5) + 1
  legend_red <- rep(1,5)
  legend_rate <- legend_value /floor(legend_value[5] )
  # legend_rate[legend_rate > 1] <- 1
  legend_blue <- 1 - legend_rate
  legend_green <- 1 - legend_rate
  par( ... )
  # add 5 red labels
  rect(legend_x1, legend_y1, legend_x2, legend_y2, col = rgb(legend_red,legend_green,legend_blue), border = NA)
  text(rep(7.5 ,5)  , seq(8,12,length.out =5) + 0.5, labels = as.character(legend_value) )
  # add maximum red labels
  rect(xleft = 2.5, ybottom = 18, xright = 5, ytop = 19, col = rgb(1,0,0), border = NA )
  text( x = 7.5 , y = 18, labels = paste("Max:\n", max(tmp_raw_matrix_vector)) )
}



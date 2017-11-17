# test tri_heatmap and rect_heatmap

tri_heatmap( as.matrix(raw_HiC_matrix[[1]][[1]][[1]][300:640,300:640 ]) )
i = 66800000/200000 - 300
j = 120200000/200000 - 300
points( i , 0, col = 'black')
points( j, 0, col = 'black')
points( x =  j/2 + i/2, y = j - i, col = 'black')
axis(1, c(0, 50, 100, 150, 200, 250, 300), labels = paste(c(0, 50, 100, 150, 200, 250, 300) /5, "M"))
title( main = 'RPMI-8226:Raw Matrix:Chr1:60M-127M')

rect_heatmap( raw_matrix = raw_HiC_matrix[[1]][[5]][[10]][1:670, 1:670] )
points( x = 34000000/200000, y = - 31800000/200000) # x:chr5 
points( x = 30000000/200000, y = - 47200000/200000) # x:chr5 , y = chr10

i = 37000000/200000 - 150
j = 120200000/200000 - 150
points( i , 0, col = 'black')
points( j, 0, col = 'black')
points( x =  j/2 + i/2, y = j - i, col = 'black')
axis(1, c(0, 50, 100, 150, 200, 250, 300), labels = paste(c(0, 50, 100, 150, 200, 250, 300) /5, "M"))
title( main = 'RPMI-8226:Raw Matrix:Chr1:60M-127M')

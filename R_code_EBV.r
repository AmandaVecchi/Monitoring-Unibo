# essential biodiversity variables = EBV

#upload image
setwd("C:/LAB/") 
library(raster) #to see the image

#raster <- #imports a single layer, yet satellite images are made of more than one layer
snt <- brick("snt_r10.tif") #also named the picture
snt #see the characteristics of the image

#plot the image
plot(snt)

#B1 = blue
#B2 = green
#B3 = red
#B4 = NIR

# R3 G2 B1
plotRGB(snt, 3, 2, 1, stretch="Lin")
plotRGB(snt, 4, 3, 2, stretch="Lin")

pairs(snt)

library(RSToolbox)
#PCA anlysis
sntPCA <- rasterPCA(snt)
sntPCA

summary(sntPCA$model) #info about output of the model
# 70% of the info
plot(sntPCA$map)

#plot RGB
plotRGB(sntPCA$map, 1, 2, 3, stretch="Lin")

###calculate standard deviation
#set the moving window
window <- matrix(1, nrow = 5, ncol = 5)

#focal --> calculates values for a neighbourhood of cells
sd_nst <- focal(sntPCA$map$PC1, w = window, fun =sd)   #we use tha map and the frst principal component; w = window; fun = function = mean, mode, max, min, sd (standar deviation)

cl <- colorRampPalette(c('dark blue','green','orange','red'))(100) 
plot(sd_snt, col=cl)

par(mfrow=c(1,2))
plotRGB(snt,4,3,2, stretch="lin")
plot(sd_snt, col=cl)















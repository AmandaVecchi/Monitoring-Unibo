## set working directory
setwd("C:/LAB/")

#importing one image
library(raster)
library(RStoolbox)
p224r63_2011 <- brick("p224r63_2011_masked.grd")

#B1 = blue
#B2 = green
#B3 = red
#B4 = NIR (leaf reflects)
#B5 = SWIR
#B6 = thermal infrared
#B7 = SWIR (it has additional properties)
#B8 = panchromatic

# RGB : 
plotRGB(p224r63_2011, r= 5, g=4, b=3, stretch="Lin")

#ggplot 
library(ggplot2)
ggRGB(p224r63_2011,5,4,3)

#SAME PROCEDURE WITH 1988 IMAGE
p224r63_1988 <- brick("p224r63_1988_masked.grd")
plotRGB(p224r63_1988, r= 5, g=4, b=3, stretch="Lin")
ggRGB(p224r63_1988,5,4,3)

#par = show both plots at the same time
par(mfrow=c(1,2))
plotRGB(p224r63_1988, r= 5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r= 5, g=4, b=3, stretch="Lin")

#see two bands together --> $ = link a part to a general part
dev.off() #remove the par
plot(p224r63_2011$B1_sre,p224r63_2011$B3_sre)

#PCA! 
#decrease resolution
p224r63_2011_res <- aggregate(p224r63_2011, fact=10)

p224r63_2011_pca <- rasterPCA(p224r63_2011_res)
#plot the map --> $ links pieces of the output
plot(p224r63_2011_pca$map) #have the default colour palette

#change colour --> colorRampPalette
cl <- colorRampPalette(c('dark grey','grey','light grey'))(100)
plot(p224r63_2011_pca$map, col=cl)
summary(p224r63_2011_pca$model)
#PC1 99.83% of the whole variation 

pairs(p224r63_2011)

#plot first three components
plotRGB(p224r63_2011_pca$map, r=1, g=2, b=3, stretch="Lin")

# same process for 1988 picture --> decrease risolution + perform PCA + plot our map 
p224r63_1988_res <- aggregate(p224r63_1988, fact=10)
p224r63_1988_pca <- rasterPCA(p224r63_1988_res) 
plot(p224r63_1988_pca$map, col=cl)
summary(p224r63_1988_pca$model)
pairs(p224r63_1988)

#difference in PCA
difpca <- p224r63_2011_pca$map - p224r63_1988_pca$map
plot(difpca)
 
cldif <- colorRampPalette(c('blue','black','yellow'))(100) 
plot(difpca$PC1,col=cldif)
 
 
 
 
 
 
 
 































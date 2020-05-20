# R_code_faPAR.r

setwd("C:/LAB/")

library(raster)
library(rasterVis)
library(rasterdiv)

plot(copNDVI)
copNDVI <- reclassify(copNDVI, cbind(253:255, NA))

levelplot(copNDVI)

#importing the file
faPAR10 <- raster("farPAR10.tif") 

levelpot(farPAR10)

#saving the plot as PDF
pfd("copNDVI.pdf")
levelplot(coNDVI)
dev.off()

pdf("faPAR.pdf")
levelplot(faPAR10)
dev.off()

#2nd day 

setwd("C:/LAB/") 
load("faPAR.RData")
library(raster)
library(rasterdiv)
library (rasterVis)

#the original faPAR from Copernicus is  2GB
# let's see how much space is needed for an 8-bit set

writeRaster( copNDVI, "copNDVI.tif")
# 5.3 MB

# to make the level plot of the faPAR that is the fraction of solar radiation absorbed by alive leaves
levelplot(faPAR10) 

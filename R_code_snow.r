setwd("~/lab/")
# setwd("/Users/utente/lab") #mac
# setwd("C:/lab/") # windows

# install.packages("ncdf4")
library(ncdf4)
library(raster)

snowmay <- raster("c_gls_SCE500_202005180000_CEURO_MODIS_V1.0.1.nc")

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 

# Exercise: plot the snow cover with the color ramp palette
plot(snowmay, c0l=cl)

# Slow manner to import the set
setwd("C:/LAB/snow/") 
snow2000 <- raster("snow2000r.tif")
snow2005 <- raster("snow2005r.tif")
snow2010 <- raster("snow2010r.tif")
snow2015 <- raster("snow2015r.tif")
snow2020 <- raster("snow2020r.tif")

par(mfrow=c(2,3))
plot(snow2000, col=cl)
plot(snow2005, col=cl)
plot(snow2010, col=cl)
plot(snow2015, col=cl)
plot(snow2020, col=cl)

# importing and plotting of many data quickly
rlist <- list.files(pattern="snow")
rlist

import <- lapply(rlist, raster)
snow.multitemp <- stack(import)
plot(snow.multitemp, col=cl)

#prediction
source("prediction.r")
plot(predicted.snow.2025.norm, col=cl)


##SECOND DAY

setwd("C:/LAB/snow/") 

#import all the snow cover images
rlist <- list.files(pattern="snow")
rlist #see the imported list

library(raster)

#lapply = applies a function over a list 
import <- lapply(rlist, raster)
snow.multitemp <- stack(import) #stack the images

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 
plot(snow.multitemp, col=cl) #plot all files

#import predicted map
  #raster function imports one image
  #brick function imports sveral bands/images all together(satellite images)
prediction <- raster("predicted.2025.norm.tif") # remember to put the file in the snow folder
plot(prediction, col = cl)

#export the output
writeRaster(prediction, "fial.tif") #writes the entire raster to a file, data  #in brackets the mname of the file

#create a PDF of the graph produced at the endo of the analyis
  #final stack
final.stack <- stack(snow.multitemp, prediction) #we took the stack of all the images (snow.multitemp), the input and stack it with the prediction we made
plot(final.stack, col=cl)

#PDF creation
pdf("Final_graph.pdf)
plot(final.stack, col=cl)
dev.off()

# PNG of the graph
png("Final_graph.png")
plot(final.stack, col=cl)
dev.off()



































































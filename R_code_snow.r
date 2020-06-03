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
setwd("~/lab/snow/")
# setwd("/Users/utente/lab/snow/") #mac
# setwd("C:/lab/snow/") # windows

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

##############

# fast version of import and plot of many data for lazy people!
rlist <- list.files(pattern="snow")
rlist

import <- lapply(rlist, raster)
snow.multitemp <- stack(import)
plot(snow.multitemp, col=cl)

# let's make a prediction
source("prediction.r")

plot(predicted.snow.2025.norm, col=cl)


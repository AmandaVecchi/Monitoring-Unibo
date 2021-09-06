insall.packages("raster")
install.packages("rasterVis")
install.packages("rasterdiv")
install.packages("rgdal")
install.packages("ggplot2")
install.packages("RStoolbox")
    
library(raster)
library(rasterVis)
library(rasterdiv)
library(rgdal)
library(ggplot2)
library(RStoolbox)

#Upload image from 2018, prefire
setwd("C:/exam/prefire/")
getwd()
rlist2018 <- list.files( pattern = "prefire")
rlist2018

#B2 = BLUE
#B3 = GREEN
#B4 = RED
#B5 = NIR
#B6 = SWIR1
#B7 = SWIR2

import_image2018 <- lapply(rlist2018, raster)
image2018 <- stack(import_image2018)
image2018

ext <- c(400000, 550000, 4350000, 4500000)
image2018_crop <- crop(image2018, ext)
plot(image2018_crop)

install.packages("raster")
install.packages("rasterVis")
install.packages("rasterdiv")
install.packages("rgdal")
install.packages("ggplot2")
install.packages("RStoolbox")
install.packages("rgeos")
install.packages("devtools")
    
library(raster)
library(rasterVis)
library(rasterdiv)
library(rgdal)
library(ggplot2)
library(RStoolbox)
library(rgeos)
library(devtools)
    
#Upload image from 2014
setwd("C:/EXAM/2014DATA/")
rlist2014 <- list.files( pattern = "2014")
rlist2014

#B2 = BLUE
#B3 = GREEN
#B4 = RED
#B5 = NIR
#B6 = SWIR1
#B7 = SWIR2

import_image2014 <- lapply(rlist2014, raster)
image2014 <- stack(import_image2014)
image2014

#crop the image
ext <- c(930000, 1050000, 4800000,4900000)
zoom(image2014, ext=ext)
image2014_crop <- crop(image2014, ext)
plot(image2014_crop)

#plot the image with RGB colours 
RGB2014 <- ggRGB(image2014_crop, r=4, g=3, b=2, stretch ="Lin")
plot(RGB2014)
    
visible2014 <- print(RGB2014 + ggtitle("Visible - Sardinia, 2014")+ theme (axis.title.x = element_blank(),axis.title.y = element_blank()))    
    
setwd("C:/project/2021/")
rlist2021 <- list.files(pattern = "2021")
rlist2021
    
import_image2021 <- lapply(rlist2021, raster)
image2021 <- stack(import_image2021)
image2021

#crop the image
ext <- c(930000, 1050000, 4800000,4900000)
zoom(image2021, ext=ext)
image2021_crop <- crop(image2021, ext)
plot(image2021_crop)

RGB2021 <- ggRGB(image2021_crop, r=4, g=3, b=2, stretch ="Lin")
plot(RGB2021) 
    
visible2021_print <- print(RGB2021 + ggtitle("Visible - Sardinia, 2021")+ theme (axis.title.x = element_blank(),axis.title.y = element_blank()))    

#plot both RGB images together
par(mfrow = c(1,2) , oma=c(0,0,2,0))
plotRBG(image2018,  r=4, g=3, b=2, stretch =”lin”)
plotRGB(image2021, r=4, g=3, b=2, stretch =”lin”)
mtext(“Visible RGB Sardinia 2018 vs Sardinia 2021”, outer = TRUE, cex = 1.5)

par(mfrow = c(1,2), oma = c(0,0,2,0))
plotRBG(image2018,  r=4, g=3, b=2, stretch =”hist”)
plotRGB(image2021, r=4, g=3, b=2, stretch =”hist”)
mtext(“Visible RGB Sardinia 2018 vs Sardinia 2021”, outer=TRUE, cex =1.5)

#CALCULATE NDVI
#NDVI = (NIR – RED) / (NIR+RED)
NDVI2018 <- (image2018[[5]] - image2018[[4]]) / (image2018[[5]] + image2018[[4]])
NDVI2018
Plot(NDVI2018, main = "NDVI 2018")

NDVI2021 <- (image2021[[5]] - image2021[[4]]) / (image2021[[5]] + image2021[[4]])
NDVI2021
Plot(NDVI2021, main = "NDVI 2021")

#PLOT BOTH NDVI IMAGES TOGETHER
par(mfrow = c(1,2), oma=c(0,0,2,0))
plot(NDVI2018, main = "NDVI 2018") 
plot(NDVI2021, main = "NDVI 2021") 
mtext(“NDVI comparison: 2018 vs 2021”, outer=TRUE, cex =1.5)

#Difference between the two years
diffNDVI <- (NDVI2018 – NDVI2021)
plot(diffNDVI)
plot(diffNDVI, main = “NDVI difference between 2018 and 2021”)

cld <- colorRampPalette(c('blue','white','red'))(100)
plot(diffNDVI, col = cld, main = "NDVI difference between 2018 and 2021")

#VISUALIZATION USING HISTOGRAMS
hist_NDVI2018 <- hist(NDVI2018, main = "Distribution of NDVI values - 2018", xlab = "NDVI", ylab = "Frequency", breaks = 30)
hist_NDVI2021 <- hist(NDVI2021, main = "Distribution of NDVI values - 2021", xlab = "NDVI", ylab = "Frequency", breaks = 30)

par(mfrow = c(1,2), oma=c(0,0,2,0))
plot(hist_NDVI2018, main = “NDVI 2018”) 
plot(hist_NDVI2021, main = “NDVI 2021”) 
mtext(“Comparison of distribution of NDVI values: 2018 vs 2021”, outer=TRUE, cex =1.5)

#calculate NBR = Normalized Burned Ratio. an index designed to highlight burnt areas 
#NBR = (NIR-SWIR) / (NIR+SWIR)
#A high NBR value indicates healthy vegetation while a low value indicates bare ground and recently burnt areas. Non-burnt areas are normally attributed to values close to zero

#2018 NBR
setwd(“C:/EXAM/2018DATA/”)
getwd()
NBR_2018 <- (image2018[[5]] - image2018[[6]]) / (image2018[[5]] + image2018[[6]])
plot(NBR_2018, main = “NBR values for 2018”)

#2021 NBR
setwd(“C:/EXAM/2021DATA/”)
getwd()
NBR_2021 <- (image2021[[5]] - image2021[[6]]) / (image2021[[5]] + image2021[[6]])
plot(NBR_2021, main = “NBR values for 2021)

#dNBR shows difference between NBR of the two years
dNBR <- (NBR_2018 – NBR_2021)
plot(dNBR, main = “Difference in NBR 2018 – 2021”)


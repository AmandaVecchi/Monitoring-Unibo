
install.packages("raster")
install.packages("rasterVis")
install.packages("rasterdiv")
install.packages("rgdal")
install.packages("ggplot2")
install.packages("RStoolbox")
install.packages("gridExtra")
    
library(raster)
library(rasterVis)
library(rasterdiv)
library(rgdal)
library(ggplot2)
library(RStoolbox)
library(gridExtra)
    
#Upload image from 2016
setwd("C:/empedrado/prefire/")
rlist2016 <- list.files( pattern = "2016")
rlist2016

#B2 = BLUE
#B3 = GREEN
#B4 = RED
#B5 = NIR
#B6 = SWIR1
#B7 = SWIR2

import_image2016 <- lapply(rlist2016, raster)
image2016 <- stack(import_image2016)
image2016

#crop the image
ext <- c(-8100000, -8000000, -4325000, -4275000)
zoom(image2016, ext=ext)
image2016_crop <- crop(image2016, ext)
plot(image2016_crop)

#plot the image with RGB colours 
RGB2016 <- ggRGB(image2016_crop, r=3, g=2, b=1, stretch ="Lin")
plot(RGB2016)
    
visible2016 <- print(RGB2016 + ggtitle("Visible - Empedrado, 2016")+ theme (axis.title.x = element_blank(),axis.title.y = element_blank()))    
    
setwd("C:/empedrado/postfire/")
rlist2019 <- list.files(pattern = "2019")
rlist2019
    
import_image2019 <- lapply(rlist2019, raster)
image2019 <- stack(import_image2019)
image2019

#crop the image
ext <- c(-8100000, -8000000, -4325000, -4275000)
zoom(image2019, ext=ext)
image2019_crop <- crop(image2019, ext)
plot(image2019_crop)

RGB2019 <- ggRGB(image2019_crop, r=4, g=3, b=2, stretch ="Lin")
plot(RGB2019) 
    
visible2019 <- print(RGB2019 + ggtitle("Visible - Empedrado, 2019")+ theme (axis.title.x = element_blank(),axis.title.y = element_blank()))    

grid.arrange(visible2016, visible2019, nrow =1)

#plot both RGB images together
par(mfrow = c(1,2) , oma=c(0,0,2,0))
plotRGB(image2016_crop, r=4, g=3, b=2, stretch ="lin")
plotRGB(image2019_crop, r=4, g=3, b=2, stretch ="lin")
mtext("Visible Empedrado 2016 vs Empedrado 2019", outer = TRUE, cex = 1.5)

par(mfrow = c(1,2), oma = c(0,0,2,0))
plotRGB(image2016_crop,  r=4, g=3, b=2, stretch ="hist")
plotRGB(image2019_crop, r=4, g=3, b=2, stretch ="hist")
mtext("Visible enhanced Empedrado 2016 vs Empedrado 2019", outer=TRUE, cex =1.5)

#CALCULATE NDVI
#NDVI = (NIR – RED) / (NIR+RED)
NDVI2016 <- (image2016_crop[[5]] - image2016_crop[[4]]) / (image2016_crop[[5]] + image2016_crop[[4]])
NDVI2016
plot(NDVI2016, main = "NDVI - Empedrado, 2016")

NDVI2019 <- (image2019_crop[[5]] - image2019_crop[[4]]) / (image2019_crop[[5]] + image2019_crop[[4]])
NDVI2019
plot(NDVI2019, main = "NDVI - Empedrado, 2019")

#PLOT BOTH NDVI IMAGES TOGETHER
par(mfrow = c(1,2), oma=c(0,0,2,0))
plot(NDVI2016, main = "NDVI 2016") 
plot(NDVI2019, main = "NDVI 2019") 
mtext("NDVI comparison: 2016 vs 2019", outer=TRUE, cex =1.5)
     
# remove cells where NDVI < 0,4
NDVI2016_mod <- reclassify(NDVI2016, cbind(-Inf, 0.4, NA))
NDVI2019_mod <- reclassify(NDVI2019, cbind(-Inf, 0.4, NA))
     
par(mfrow=c(1,2))
plot(NDVI2016_mod, main="Vegetation 2016",  axes=FALSE)
plot(NDVI2019_mod, main="Vegetation 2019", axes=FALSE)

#Difference between the two years
diffNDVI <- (NDVI2016 - NDVI2019)
plot(diffNDVI)
plot(diffNDVI, main = "NDVI difference between 2016 and 2019")

cld <- colorRampPalette(c('blue','white','red'))(100)
plot(diffNDVI, col = cld, main = "NDVI difference between 2016 and 2019")

#VISUALIZATION USING HISTOGRAMS
hist_NDVI2016 <- hist(NDVI2016, main = "Distribution of NDVI values - 2016", xlab = "NDVI", ylab = "Frequency", breaks = 50)
hist_NDVI2016 <- hist(NDVI2016, main = "Distribution of NDVI values - 2016", xlab = "NDVI", ylab = "Frequency", breaks = 50, col = "light pink",border = "black")
     
hist_NDVI2019 <- hist(NDVI2019, main = "Distribution of NDVI values - 2019", xlab = "NDVI", ylab = "Frequency", breaks = 50)
hist_NDVI2019 <- hist(NDVI2019, main = "Distribution of NDVI values - 2019", xlab = "NDVI", ylab = "Frequency", breaks = 50, col = "light blue",border = "black")
     
par(mfrow = c(1,2), oma=c(0,0,2,0))
plot(hist_NDVI2016, main = "NDVI 2016", col = "light pink",border = "black") 
plot(hist_NDVI2019, main = "NDVI 2019", col = "light blue",border = "black") 
mtext("Comparison of distribution of NDVI values: 2016 vs 2019", outer=TRUE, cex =1.5)
     
col2rgb("lightblue")
col2rgb(c("lightblue", "lightgreen", "pink"))
mycol <- rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50") #make the color transparent by 50%
col2rgb("lightblue") #to see the red, green and blue values you need 
    ## red    173
    ## green  216
    ## blue   230
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink") 

par(mai=rep(0.5, 4)) #set the margins for the image
layout(matrix(c(1,1,2,2,0,3,3,0), ncol = 4, byrow = TRUE)) #devide the plotting space
plot(hist_NDVI2016, col=c2, main="NDVI 2016", xlab = "NDVI")
plot(hist_NDVI2019, col=c1, main="NDVI 2019", xlab = "NDVI")
plot(hist_NDVI2016, col = c2, xlim = c(-0.5, 1), main="Comparison between NDVI",xlab = "NDVI" )
plot(hist_NDVI2019, add = TRUE, col = c1)
     
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

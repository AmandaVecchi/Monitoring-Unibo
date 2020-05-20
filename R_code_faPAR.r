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

###regression model between faPAR and NDVI
erosion <- c(12, 14, 16, 24, 26, 40, 55, 67) # amount of erosion in a certain area
hm <- c(30, 100, 150, 200, 260, 340, 460, 600) #heavy metals in the area, nichel for example. In ppm

#plot between these two variables
plot(erosion, hm, col = "red", pch = 19, xlab = "Erosion", ylab = "Heavy Metals")
plot(erosion, hm, col = "red", pch = 19, xlab = "Erosion", ylab = "Heavy Metals", cex = 2) #exagerate dimension of poin

#linear model between the two variables
model1 <- lm(hm ~ erosion) #first we named the model
summary(model1) #gives info about the model --> R^2 ranges between 1 and 0, the higher the better the variables are significantly related = pattern we observe is not random. p-valus is small, the pattern is not random. 
#variables are related

#plot the linear regression, the straight line
abline(model1)

#faPAR vs NDVI model
#load the faPAR data previously created
#set working directory
setwd("C:/LAB/")
library(raster) #to recognize the picture of faPAR
library(rasterdiv)
faPAR10 <- raster("faPAR10.tif") 
plot(faPAR10)
plot(copNDVI)

#remove values, water data
copNDVI <- reclassify(copNDVI, cbind(253:255, NA), right=TRUE)

#two variables: copNDVI and faPAR. see how much they are related --> new model

#see number of cells in the raster = ncells
faPAR10 

#function defined to directly select random points from an image. x=raster file, n=number of points we want to select
install.packages("sf")
library(sf) # to call st_* functions
random.points <- function(x,n)
{
lin <- rasterToContour(is.na(x))
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') # st_union to dissolve geometries
pts <- spsample(pol[1,], n, type = 'random')
}

pts <- random.points(faPAR10, 10000) #firstly we named the function, then create it

#extract values fo faPAR and NDVI
copNDVIp <- extract(copNDVI, pts) #extract data and put them on a point. for each point we have the value of the variable
faPAR10 <- extract(faPAR10, pts)

copNDVIp #see the vlues of the points chose by random procedure. NA = points in the sea
faPARp #see values for eac point

#photosynthesis vs. biomass
model2 <- lm(faPAR10 ~ copNDVI)
summary(model2)

#regression model for the relationship between the two variables
plot(copNDVIp, faPARp, col = "green", xlab = "Biomass", ylab = "Photosynthesis")
abline(model2, col = "red")
















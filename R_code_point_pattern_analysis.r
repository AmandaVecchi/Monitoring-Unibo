## Point pattern analysis: density maps

install.packages("spatstat")
library(spatstat)

#attach our dataset, the one we are going to use
attach(covid)

#what are the coordinates of the datases and what are the extensions of these cooridnates
head(covid)
covids <- ppp(lon, lat, c(-180,180), c(-90,90))  #ppp= planar point pattern, specify the variables and the range

#density map
d <- density(covids)

#show density map
plot(d)

#put points on the density map
points(covids)

setwd("C:/LAB/")
load(".RData")
ls()

#covids: point pattern
#d = density map
library(spatstat)

plot(d)
points(covids)

install.packages("rgdal")
library(rgdal)

# let's input vector ines (x0y0, x1y1, x2y2..) #download coastline file from IOL and put it in older LAB
coastline <- readOGR("ne_10m_coastline.shp")

plot(coastline, add=T) #adds coastline to previous plot of covids

# change of the colour, we are making ramp palette   c=array of colours   (100): all the possible colours from yellow to red
cl <- colorRampPalette(c("yellow","orange","red"))(100)
plot(d, col=cl)
points(covids)
plot(coastlines, add=T)

# Exercise: new colour ramp palette
clr <- colorRampPalette(c("light green","yellow","orange","violet"))(100)
plot(d, col=clr, main="densitites of covid-19")
points(covids)
plot(coastlines, add=T)

# export as pdf or png("covid_density.png")
pdf("covid_density.pdf")
clr <- colorRampPalette(c("light green","yellow","orange","violet"))(100)
plot(d, col=clr, main="densitites of covid-19")
points(covids)
plot(coastlines, add=T)
dev.off()











# R code to view biomass over the world and calculate changes in ecosystems functions
# energy
# chemical cycling
# proxies

install.packages("rasterdiv")
install.packages("rasterVis") #for visualization

library(rasterdiv)
library(raserVis)

data(copNDVI) #Declare dataset --> NDVI = NIR + RED
plot(copNDVI)

#reclassifying data --> create a new set
copNDVI <- reclassify(copNDVI, cbind(253:255, NA))
levelplot(copNDVI) # Mean biomass over last 30 years, alle eco functions and services that plants are giving

#aggregate values --> change scale, resolution 
copNDVI10 <- aggregate(copNDVI, fact=10) #increase pixel dimension and give a name to the image
levelplot(copNDVI10)

copNDVI100 <- aggregate(copNDVI, fact=100)
levelplot(copNDVI100)

## uploading images for deforestation
#set working directory
setwd("C:/LAB/")
library(raster)
defor1 <- brick("defor1_.jpg")
defor2 <- brick("defor2_.jpg")

# plot two images --> 
band 1= NIR
band2 = RED
band3 = GREEN

plotRGB(defor1, r=1, g=2, b=3, stretch="Lin") 
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")

#lot together
par(mfrow=c(1,2))
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")

#calculate the amount of loss of DVI
dvi1 <- defor1$defor1_.1 - defor1$defor1_.2
dvi2 <- defor2$defor2_.1 - defor2$defor2_.2

#plot together
cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100) # specifying a color scheme
par(mfrow=c(1,2))
plot(dvi1, col = cl)
plot(dvi2, col = cl)

difdvi <- dvi1 - dvi2

dev.ogg() #to remvove previous part
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(difdvi, col = cld)

#histogram showing freq of the maps
hist(difdvi)




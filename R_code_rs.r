#R code for remote sensing data and analysis

#Setting the working directory, uploading packages we are using and state we are using them
setwd("C:/LAB/")  
install.packages("raster")
install.packages("RStoolbox")
#RStoolbox is a package providing tools for remote sensing processing. 
library(raster) 
library(RStoolbox)

#The "brick" function allows us to import images, files
p224r63_2011 <- brick("p224r63_2011_masked.grd") #Here we also give a name to the imported image
plot(p224r63_2011)

#The "colorRampPalette" is a function that takes a set of given colours and returns a character vector of colours
#Creating a new colorRampPalette anche also naming it
cl <- colorRampPalette(c('red','blue'))(100) 
plot(p224r63_2011, col=cl) #Plotting the image using the newely created palette


###LANDSAT is a constellation of satellites orbiting around the Earth and collecting data used to study different environments and their changes in time
# Bands of landsat
# B1: blue band
# B2: green band
# B3: red band
# B4: NIR (infrared) band 
#These bands can be plotted in different ramp palette

#In order to close the previous plot and continue working
dev.off()
 
#The "par" function creates a multiframe of different plots in the dataset
#The "mfrow" argument is used to set the parameters 
par(mfrow=c(2,2)) #We set the number of columns - 2- and rows - 2.

#Creating a plot for the band B1 using a specific palette
clb <- colorRampPalette(c('dark blue','blue','light blue'))(100) 
plot(p224r63_2011$B1_sre, col=clb)   #The $ symbol is used to link every single band to an image

#Creating a plot for the band B2 using a specific palette
clg <- colorRampPalette(c('dark green','green','light green'))(100) 
plot(p224r63_2011$B2_sre, col=clg)

#creating a plot for band B3 using a specific palette
clr <- colorRampPalette(c('dark red','red','pink'))(100) 
plot(p224r63_2011$B3_sre, col=clr)

#Creating a plot for band B4 using a specific palette
cln <- colorRampPalette(c('red','orange','yellow'))(100) 
plot(p224r63_2011$B4_sre, col=cln)

#Creating a graph showing all 4 images created before 
par(mfrow=c(4,1))

###RGB plotting
#In the computer there are 3 components that make the colors visible: RGB system
#R --> we associate the 3rd red band to R
#G --> we associate the 2nd green band to G
#B --> we associate the 1st blu band to B

dev.off() #To close the previous graph

#The "plotRGB" function creates a Red-Green-Blue plot based on three layers
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")
#In the argument state the image we want to plot and also the correspondence between the RGB system and the bands of landsat
#stretch = "Lin" is used to stretch the values to increase the contrast of the image

# Exercise: NIR on top of the G component of the RGB 
# invert 4 with 3
plotRGB(p224r63_2011, r=3, g=4, b=2, stretch="Lin") 
plotRGB(p224r63_2011, r=3, g=2, b=4, stretch="Lin")


###Working on the 1988 image

#Setting the working directory, loading the previously used workspace and see the packages and function used
setwd("C:/LAB/")
load("First.RData")
ls()

library(raster)
p224r63_1988 <- brick("p224r63_1988_masked.grd") #Importing the 1988 image and naming it
 
#Plotting with the RBG system both images: one from 2011 and the other from 1988
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=3, g=2, b=1, stretch="Lin") 
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")

#Plotting the images in false colours RGB - 432
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="Lin") 
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")

#Enhancing the noise of the image
par(mfrow=c(2,1))
plotRGB(p224r63_1988, r=4, g=3, b=2, stretch="hist") 
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="hist")

#DVI for the two years: compare with a difference in time
#DVI = NIR - RED = difference vegetation index. This index is based on how differently plants reflect light. Distinguishes between soil and vegetation
#NDVI = (NIR - RED) / (NIR + RED); Satellite indicator of the presence of vegetation on the Earth's surface and its changes in time.

#DVI for 2011
dvi2011 <- p224r63_2011$B4_sre - p224r63_2011$B3_sre 
cl <- colorRampPalette(c("yellow","light blue","lightpink4"))(100) 
plot(dvi2011, col=cl)

#DVI for 1988
dvi1988<- p224r63_1988$B4_sre - p224r63_1988$B3_sre 
cl <- colorRampPalette(c("yellow","light blue","lightpink4"))(100) 
plot(dvi1988, col=cl)

#Plotting both images together
par(mfrow=c(2,1))
plot(dvi1988)
plot(dvi2011)

#Plotting both graphs together with a new palette
par(mfrow=c(2,1))
cldvi <- colorRampPalette(c('red','orange','yellow'))(100)  
plot(dvi1988, col=cldvi)
plot(dvi2011, col=cldvi)

#Plotting the difference from one year to the other to see what changed
diff <- dvi2011 - dvi1988 #naming the new image
plot(diff)

#Plotting the difference between 2011 and 1988 with a new palette
cldif <- colorRampPalette(c('blue','white','red'))(100) #
plot(diff, col=cldif)

#Changing the grain of our image, the size of the pixels
#The "aggregate" function changes the dimension of the pixels of the required units
#Smaller pixels can help us get better info
p224r63_2011res <- aggregate(p224r63_2011, fact= 10) #fact= 10 states the factor of aggregation of pixels
p224r63_2011res100 <- aggregate(p224r63_2011, fact= 100) 

#plot both together
par(mfrow=c(3,1))
plotRGB(p224r63_2011, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res100, r=4, g=3, b=2, stretch="Lin")
 






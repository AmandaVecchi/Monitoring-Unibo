setwd("C:/LAB/NO2/")

#to see image
library(raster)

##import all the NO2 data by lapply function
#create a list of objects that have the same pattern (part of the name)
rlist <- list.files(pattern = "EN")
rlist
               
#lapply applies a function over a list of object
import <- lapply(rlist, raster)

#stack, creates a stack of all images
EN <- stack(import)

#plot
cl <- colorRampPalette(c('red','orange','yellow'))(100) 
plot(EN)

# dev.off() #if you have trouble visualizing the following plot

#plot only 2 of the images (first and last one)
par(mfrow = c(1,2))
plot(EN$EN_0001, col = cl)
plot(EN$EN_0013, col = cl)

dev.off()

# RGB space
plotRGB(EN, r=1, g=7, b=13, stretch = "lin") #each colour is associated with the number of the image in the stack
#where we see red the first image has high values and so on

#difference map
dif <- EN$EN_0013 - EN$EN_0001
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(dif, col=cld)

#BOXplot = produces a box-and-whiskers plot
#quantitative estimate
boxplot(EN) #of the whole stack

#remove outliers
boxplot(EN, outline=F)
 
#plot horizontally
boxplot(EN, outline=F, horizontal=T)

# dev.off()
#name the axes
boxplot(EN, outline=F, horizontal=T, axes=T) #median is quite similar in all boxplot, but max value is decreasing
 
#if there is a decrease in NO2,diff values should lay below the y=x line
  #plot all the images, all the pixels
plot(EN$EN_0001, EN$EN_0013) 
  #add the 1:1 line
plot(EN$EN_0001, EN$EN_0013)
abline(0,1,col="red")

setwd("C:/LAB/snow/")
rlist <- list.files(pattern="snow20")
rlist

import <- lapply(rlist, raster)
snow.multitemp <- stack(import)
snow.multipemp

plot(snow.multitemp$snow2010r, snow.multitemp$snow2020r)
abline(0,1, col="red")

plot(snow.multitemp$snow2000r, snow.multitemp$snow2020r)
abline(0,1, col="red")











 





































# R code for spatial view of points

# show the points in a map
library(sp)

data(meuse)

#Look at the data. part of the data, the first 6 lines
head(meuse)

# this is a spatial dataset. Tell R this dataset contains coordinates
coordinates(meuse) = ~ x+y    # = ~ specify which are the coordinates, firstly we have to specify the dataset that contains the coordinates

plot(meuse) # position of coordinates

# plot the spatial amount of zinc
spplot(meuse, "zinc") # declare the set and declare the variable

# plot the spatial amount of copper 
spplot(meuse, "copper", main="copper concentration") #change title --> main="..."

# show the coordinates in bubbles with different size depending on the concentration in that point
bubble(meuse, "zinc")
bubble(meuse, "zinc", main= "Zinc concentration")

#bubble copper in red
bubble(meuse, "copper", col="red")
bubble(meuse, "copper", main= "Copper concentration", col= "red")

# take data from another source and upload them into R. --> download covid.agg
#Put covid file into the folder LAB in our computer

#setting the working directory: LAB
#windows: 
setwd("C:/LAB/")
mac:
setwd("/users/nameofthecomputer/LAB")

#give a name to the dataset = covid. link the object (our dataset) to function.
covid <- read.table("covid_agg.csv", head=TRUE) # use quotes because it is a table outside of R. head=TRUE or head=T to tell R that the first line is the header of the table
head(covid)

#spatial plot --> number of cases per country
#attach dataset to R 
attach(covid)
plot(country,cases)

#if dataset is not attach --> plot(covid$country,covid$cases)

#las = how to direction lables of the axis
plot(country, cases, las=0) #parallel lables
plot(country, cases, las=1) #horizontal lables
plot(country, cases, las=2) #perpendicular lables
plot(country, cases, las=3) #vertical lables

#reduce sixe of characters of axis  
plot(country, cases, las=3, cex.axis=0.5)
plot(country, cases, las=3, cex.axis=0.7)

#install ggplot2 package
install.packages("ggplot2")
library(ggplot2)

#load previously saver R data
#setting woring directory
setwd("C/LAB/)
load(".RData") # put the name of the workspace saved
ls() 

#use ggplot2 package
library(ggplot2)

data(mpg)
head(mpg) 
#key components og ggplot2: data, aestetics, geometry #aes = variables, geometry is separated
ggplot(mpg, aes(x=displ,y=hwy)) + geom_point()
ggplot(mpg, aes(x=displ,y=hwy)) + geom_line() #change the geometry, to bars, lines 
ggplot(mpg, aes(x=displ,y=hwy)) + geom_polygon()

#ggplot for our covid dataset
head(covid)
ggplot(covid, aes(x=lon, y=lat, size=cases)) + geom_point() #size =dimension of point changes with respect to a variable





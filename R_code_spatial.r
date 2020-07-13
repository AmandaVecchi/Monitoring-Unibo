## R code for spatial view of points

#Use the spatial library in order to show the points in a map
library(sp)

#Declare the database we are going to use
data(meuse)

#head = shows the first six lines of our database
head(meuse)

#coordinates = sets spatial coordinates. We are declaring R this is a spatial dataset that contains coordinates
# ~ specifies which are the coordinates
coordinates(meuse) = ~ x+y   

#plot = function for plotting objects
plot(meuse)

#ssplot = plots several layers with a single legend for all maps
#In this case we are plotting the spatial amount of zinc
spplot(meuse, "zinc") #declare the set and declare the variable

#Plot the spatial amount of copper 
#main = titles the plot
spplot(meuse, "copper", main="copper concentration")

#bubble = creates a bubble plot of spatial data
#Each bubble is a coordinate and the size depends on the concentration in that point. Bigger the bubble, higher the concentration
bubble(meuse, "zinc")
bubble(meuse, "zinc", main= "Zinc concentration")

#Pplot the copper concentration with a bubble plot in which the bubbles are colored in red
bubble(meuse, "copper", col="red")
bubble(meuse, "copper", main= "Copper concentration", col= "red")

## Analysis of data regarding COVID-19
#Download a new database called "covid.agg" and put the file into the folder "LAB" in our computer

#Setting the working directory: LAB
setwd("C:/LAB/")

#read.table = reads a file in table format and creates a data frame from it
#arrow = give a name to the dataset
#quotes are used for each object that is outside of R and in getting imported.
covid <- read.table("covid_agg.csv", head=TRUE)
# head=TRUE or head=T tells R that the first line is the header of the table

#Show the first lines of our new dataset
head(covid)

#attach = attaches the databae to the search path. Objects in the database can be accessed by simply giving their names
attach(covid)
plot(country,cases) #plot the number of cases per country

##If the dataset is not attaced we use $ to refer to a specific object relative to a specific data frame
##plot(covid$country, covid$cases)

#las = numeric value indicating the orientation of the labels and any other text added to a plot 
plot(country, cases, las=0) #lables parallel to axis
plot(country, cases, las=1) #lables horizontal to axis
plot(country, cases, las=2) #lables perpendicular to axis
plot(country, cases, las=3) #lables vertical to axis

#cex = number indicating the amount by which plotting text and symbols should be scaled relative to the default
#1 = default
#1.5 is 50% larger
#0.5 is 50% smaller
#cex.axis	= changes the size of the axis annotation
plot(country, cases, las=3, cex.axis=0.5) #reducing size of 50%
plot(country, cases, las=3, cex.axis=0.7) #reducing size of 70%

########

#Install "ggplot2" package
#ggplot2 is a system for creating graphics
install.packages("ggplot2")
library(ggplot2) #in order to make use of the package

#setting woring directory
#Load previously saved R data
setwd("C:/LAB/)
load("First.RData") 
#First.RData is the name of the saved working space

#ls = shows what data sets and functions we defined and used in the recalled workingspace
ls() 
library(ggplot2)

data(mpg)
head(mpg) 

#key components of ggplot2: data, aesthetics, geometry 
#aes = declare the variables
#geom_point() = declare the geometry of the plot 
ggplot(mpg, aes(x=displ,y=hwy)) + geom_point() 
ggplot(mpg, aes(x=displ,y=hwy)) + geom_line() #change the geometry to lines 
ggplot(mpg, aes(x=displ,y=hwy)) + geom_polygon() #change the geometry to polygons

#ggplot for our covid dataset
head(covid)
ggplot(covid, aes(x=lon, y=lat, size=cases)) + geom_point() 
#size = x -->  dimension of points changes with respect to a variable x
#In this case the size of point changes accordingly to the number of cases





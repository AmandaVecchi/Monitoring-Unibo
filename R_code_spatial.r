## R code for spatial view of points

#Recalling the previously intalled package to import, manipulate and export spatial data
library(sp)
data(meuse)
head(meuse)

#The "coordinates" function estimates spatial coordinates for the elements in the dataset.
#We are declaring R this is a spatial dataset that contains points related to specific coordinates
coordinates(meuse) = ~ x+y #The tilde operator specifies which are the coordinates
plot(meuse)

#The "spplot" function allows us to plot several layers with a single legend for all maps
#In this case we are plotting the spatial amount of zinc
spplot(meuse, "zinc") #In the argument we declare the dataset and the variable we are plotting
spplot(meuse, "copper", main="copper concentration") #Main allows us to title the plot we create

#The "bubble" function creates a bubble plot of spatial data
#Each bubble represent a coordinate and a variable; the size of the bubble depends on the concentration in that point. Bigger the bubble, higher the concentration
bubble(meuse, "zinc")
bubble(meuse, "zinc", main= "Zinc concentration")

bubble(meuse, "copper", col="red")
bubble(meuse, "copper", main= "Copper concentration", col= "red")


## Analysis of data regarding COVID-19
#Download a new database called "covid.agg" and put the file into the folder "LAB" in our computer

#Setting the working directory (where we put the images or tables we are using for our analysis): LAB
setwd("C:/LAB/")

#The "read.table" function reads a file in table format and creates a data frame from it. Rows correspond to cases and columns to variables
#<- is an operator that assigns a name to the function we perform
covid <- read.table("covid_agg.csv", head=TRUE) 
#quotes are used for each object that is outside of R and in getting imported.
#head=TRUE or head=T tells R that the first line in the table is the header of the table

head(covid)
attach(covid)
plot(country,cases) #We are plotting the number of cases per country

##If the dataset is not attaced we use $ to refer to a specific object relative to a specific data frame
##plot(covid$country, covid$cases)

#las is a numeric value indicating the orientation of the labels and any other text added to a plot 
plot(country, cases, las=0) #lables parallel to axis
plot(country, cases, las=1) #lables horizontal to axis
plot(country, cases, las=2) #lables perpendicular to axis
plot(country, cases, las=3) #lables vertical to axis
plot(country, cases, las=3, cex.axis=0.5) #reducing size of 50%
plot(country, cases, las=3, cex.axis=0.7) #reducing size of 70%

########

#Install and stating we are making use of the "ggplot2" package
#The ggplot2 package allows us to plot variables 
install.packages("ggplot2")
library(ggplot2) 

#Setting woring directory: LAB
#The "load" function allows us to upload previously saved R datasets
setwd("C:/LAB/")
load("First.RData") #In brackets there is the name of the saved working space

#The "ls" function shows what datasets and functions we defined and used in the recalled working space
ls() 
data(mpg) #We are using a new dataset
head(mpg) 

#Key components of ggplot2: data, aesthetics, geometry 
#aes = declare the variables
#geom_point() = declare the geometry of the plot 
ggplot(mpg, aes(x=displ,y=hwy)) + geom_point() 
ggplot(mpg, aes(x=displ,y=hwy)) + geom_line() #change the geometry to lines 
ggplot(mpg, aes(x=displ,y=hwy)) + geom_polygon() #change the geometry to polygons

#Let's create a ggplot for our covid dataset
head(covid)
ggplot(covid, aes(x=lon, y=lat, size=cases)) + geom_point() 
#size = cases -->  dimension of points changes with respect to the variable we stated
#In this case the size of points changes accordingly to the number of cases





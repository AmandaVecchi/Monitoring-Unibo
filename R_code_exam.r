#R_code_exam.r

# 1. R code first
# 2. R code spatial 
# 3. R code multipanel 
# 4. R code point pattern analysis
# 5. R code multivariate analysis 
# 6. R code remote sensing
# 7. R code ecosystem function
# 8. R code pca remote sensing 
# 9. R code radiance
# 10. R code faPAR
# 11. R code EBV
# 12. R code NO2 
# 13. R code snow
# 14. R code crop
# 15. R code interpolation
# 16. R code sdm
# 17. R code exam project 


########## 1. R code first

install.packages("sp") # shows the points in a map
library(sp) # attach an add-on package
data(meuse) # to load specific dataset
meuse # see how the meuse dataset is structured
head(meuse) # shows the first rown in the dataset

#Plotting to variables to see if their concentration is related (zinc - copper)

attach(meuse) #the attach function allows to use the variables of a data.frame without always recalling the dataset
plot(zinc,coper) # plotting zinc vs copper concentration
plot(zinc,copper,col="green") # "col" allows to modify the color of the plotted symbols. IMPO: colours in brackets
plot(zinc,copper,col="green",pch=19) # "pch" is used to specify the type of point we use
plot(zinc,copper,col="green",pch=19,cex=2) # "cex" controls the size of the symbol: 1=default, 1.5 = 50%largers, 0.5 = 50% smaller

############################################################
############################################################
############################################################

########## 2. R code spatial

install.packages("sp")
library(sp)
data(use)
head(meuse)

coordinates(meuse)= ~ x+y # "coordinates" convert values in "meuse" into spatial coordinates equal the group of X and Y
plot(meuse) # plotting the points as coordinates
spplot(meuse, "zinc")  #declare which dataset we want to use and add a variable 

#EXERCISE: plot the spatial amount of copper
spplot(meuse, "copper")
spplot(meuse, "copper", main="Copper Concentration") # "main" adds a title to the plot. IMPO: the title has to be in brackets
bubble(meuse, "zinc", main="Zinc Concentration") # "bubble" plots the points with bubbles whose size reflects its value

#EXERCISE: bubble copper in red
bubble(meuse, "copper", col="red", main="Copper Concentration")

###### NEW EXERCISE. IMPORTING NEW DATA
# download the covid_agg.csv file into a folder called "LAB" into C:

setwd("C:/LAB/") # setting the working directory = place where we keep all the files and elements we use for a project
install.packages("ggplot2") #new package to plot spatially 
library(ggplot2)

# importing the downloaded table into R as a data.frame
covid <- read.table("covid_agg.csv", head=TRUE ) # "head = TRUE" or "head = T" states that the first is the title, not a line of data. 
head(covid)
attach(covid)
plot(country,cases) #plotting the number of cases for each country

#plot(covid$country, covid$cases) in this case we specify the variables "attaching" them to the dataset they belong to

#this graph doesn't show all the countries' name, so we have to change the direction of the lables
plot(country,cases, las=0) #las=0 parallel labels to axes
plot(country,cases, las=1) #las=1 horizontal labels to axes
plot(country,cases, las=2) #las=2 perpendicular labels to axes
plot(country,cases, las=3) #las=3 vertical labels to axes
plot(country,cases, las=3, cex.axis=0.5) #cex.axis changes the size of the axes label with

data(mpg)
head(mpg)

###### let's continue

load("spatial.RData") #upload in R the previously saved project
ls() #list shows all the data and variables we have used so far in the project

ggplot(mpg, aes(x=displ,y=hwy))+geom_point() #(dataset we are using, aes(x=name of the variable we plot, y=name of the variable we plot) ) + : geometry we want to use to plot that data
ggplot(mpg, aes(x=displ,y=hwy))+geom_line() # "geom_line creates a linear graph
ggplot(mpg, aes(x=displ,y=hwy))+geom_polygon() # "geom_polygon plots data with polygonal geometry

###### ggplot of covid data
head(covid)
ggplot(covid, aes(x=lon,y=lat,size=cases))+geom_point()  # "size"= the size of the point changes according to the number of cases in each country

############################################################
############################################################
############################################################

########## 3. R code multipanel

install.packages("sp")
install.packages("GGally") # "GGally" is an extension of ggplot2  
library(sp) #require(sp) is the same 
library(GGally) 

data(meuse) 
attach (meuse)

#EXERCISE: see the names of the variabes and plot cadmium versus zinc
head(meuse) # names (meuse) will do the same job
plot(cadmium,zinc, pch=15, col="red", cex=2) #points plotted in red and size is magnified x2

#EXERCISE: make all the possible pairwise plots of the dataset
#plot(x,cadmium)
#plot(x,zinc)
#plot...
#plot isn't a god idea; we use another function: "pairs" that creates a matrix of scatterplots with all possible pairwise plots
pairs(meuse) #the result is a multipanel showing multiple graphs in one figure

pairs(~ cadmium+copper+lead+zinc, data=meuse)  #grouping variabiles
# "~" is called "tilde" -> blocnum alt + 126 
# pairs(meuse[,3:6]) #name of dataset + make the subset of the dataset  "," comma means "start from" ":" means "until"

#EXERCISE: prettify this graph 
pairs(meuse[,3:6], col="blue", pch=4, cex=1)
ggpairs(meuse[,3:6])

############################################################
############################################################
############################################################

########## 4. Rcode point pattern analyses (density maps)

#Installing the package for statystical analyses of spatial point pattern = Spatial Point Pattern Analysis
install.packages("spatstat") 
library(spatstat)

attach(covid)
head(covid)

#naming the object we are going to prepara
covids <-ppp(lon,lat,c(-180,180), c(-90,90))  # "ppp" creates a point pattern dataset in a 2D plane
#ppp(variable x, variable y, c(range of the x variable), c(range of the y variable)

### Density function
d <- density(covids)
plot(d) #sowing the density map
points(covid)

##### let's continue
setwd("C:/LAB")
load("pointpattern.RData")
ls()

install.packages("rgdal")
library(rgdald)
library(spatstat) 

plot(d)
points(covids)

#Plot the coastlines on the graph to have a better representation of cases in the world
#Importing the file containing the coastline info
coastlines<- readOGR("ne_10m_coastline.shp") # "readOGR" transforms a OGR data source and layer into a suitable Spatial vector object

plot(d)
points(covids)
plot(coastline, add=T)

#Let's prettyfy the graph
#create a new palette of colours with a specific range of colours
cl <- colorRampPalette (c("yellow","orange","red"))(100) # 100=shades within the range of colours

plot(d, col=clr, main="Densities of covid-19") #Adding a title to the map
plot(d,col=cl)
points(covids)
plot(coastlines, add=T)

# EXERCISE: create a new colour ramp palette and plot again the graph
cl <- colorRampPalette (c("blue","green","yellow"))(100)
plot(d,col=cl)
points (covids)
plot(coastlines, add=T)

#Export the graph we created in -PDF
pdf("covid_density.pdf")  #or png(covid_density-png)

#Close the project
dev.off () #dev=device 

############################################################
############################################################
############################################################

########## 5. R code for multivariate analyses

install.packages("vegan") # "vegan" is a community ecology package
library(vegan)
setwd("C:/LAB/")

#import the table 
biomes <- read.table("biomes.csv", head=T, sep=",") #in the biomes there is an header (so T=true), the values are separated by comma
head(biomes) #or view(biomes) -> to view the dataset 

#Import a table 
biomes <- read.table("biomes.csv", head=T, sep=",") 
#In the biomes table different columns are separated by commas, we have to tell R this info -> sep= ","
head(biomes) #or view(biomes) 

#Let's see how species are related to eachother
#use the decorana function = DEtrended CORrespondence ANAlysis
multivar <- decorana(biomes)
plot(multivar)

#Import an additional useful table
biomes_types <- read.table("biomes_types.csv", head=T, sep=",")
head(biomes_types)
attach(biomes_types) #declaring we are working with this new data.frame













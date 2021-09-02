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

#Drawing an ellipse that connect all the points belonging to the same biome
ordiellipse(multivar, type, col=1:4, kind = "ehull", lwd=3)
# kind= type of graph --> hull is a convex shape and "e" is for ellipse
# 4 different biomes, so 4 different colors or we can write col=c("green","blue","red","black")

ordispider(multivar, type, col=1:4, label=T) #create a spider diagram

############################################################
############################################################
############################################################

########## 6. R code for remote sensing (RS)

setwd("C:/LAB/") 
install.packages("raster") # raster allows us to work with images
install.packages("RStoolbox") # package for remote sensing image, processing and analysis such as calculating spectral indeces, principal component transformation
library(raster) 

#Let's import an image
p224r63_2011 <- brick("p224r63_2011_masked.grd") 
plot(p224r63_2011)

#Plot the image with different colours -> create e new palette
cl <- colorRampPalette(c('black','grey','light grey'))(100) 
plot(p224r63_2011, col=cl)

#EXERCISE: plot the image with a new color ramp palette
cl <- colorRampPalette(c('red','blue'))(100) 
plot(p224r63_2011, col=cl)

# Bands of landsat
# B1: blue band
# B2: green band
# B3: red band
# B4: NIR (infrared) band (just after the red in the electromagnetic spectrum)

#Closing the graph
dev.off()

#Create an image with different plots
par(mfrow=c(2,2)) #number of row x column
	
# B1: blue band
clb <- colorRampPalette(c('dark blue','blue','light blue'))(100) 
plot(p224r63_2011$B1_sre, col=clb)     #$ is used to link every single band to an image
	
# B2: green band
clg <- colorRampPalette(c('dark green','green','light green'))(100) 
plot(p224r63_2011$B2_sre, col=clg)  
	
# B3: red band
clr <- colorRampPalette(c('dark red','red','pink'))(100) 
plot(p224r63_2011$B3_sre, col=clr)  
	
# B4: NIR band
cln <- colorRampPalette(c('red','orange','yellow'))(100) 
plot(p224r63_2011$B4_sre, col=cln)  

#Create an image with all the images above arranged in 1 column
par(mfrow=c(4,1))
dev.off()

# in the computer there are 3 components that make the colors visible: RGB system
# R --> we associate the 3rd red band to R
# G --> we associate the 2nd green band to G
# B --> we associate the 1st blu band to B

#plotRGB
plotRGB(p224r63_2011, r=3, g=2, b=1, stretch="Lin")
#plotRGB(name of the images, correspondence between the RGB system and the bands of landsat, the color is stretched linearily 

# EXERCISE: NIR on top of the G component of the RGB 
# invert 4 with 3
plotRGB(p224r63_2011, r=3, g=4, b=2, stretch="Lin") 
plotRGB(p224r63_2011, r=3, g=2, b=4, stretch="Lin")

###### Let's continue
setwd("C:/LAB/")
load("rs.RData")
ls()
library(raster)

#Working with other images, older images
p224r63_1988_masked <- brick("p224r63_1988_masked.grd")
plot(p224r63_1988_masked)

# EXERCISE: plot in visible RGB 321 both images (2011 and 1988)
par(mfrow=c(2,1))
plotRGB(p224r63_1988_masked, r=3, g=2, b=1, stretch="Lin")
plotRGB(p224r63_2011_masked, r=3, g=2, b=1, stretch="Lin")

# EXERCISE: plot in visible RGB 432 both images
par(mfrow=c(2,1))
plotRGB(p224r63_1988_masked, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011_masked, r=4, g=3, b=2, stretch="Lin")

#let's enhance the noise of the images
par(mfrow=c(2,1))
plotRGB(p224r63_1988_masked, r=4, g=3, b=2, stretch="hist")
plotRGB(p224r63_2011_masked, r=4, g=3, b=2, stretch="hist")

#DVI= NIR - RED --> stressed plants have very low value of difference vegetation index
dvi2011 <- p224r63_2011_masked$B4_sre - p224r63_2011_masked$B3_sre #calculate DVI -> NIR is B4 - B3 in landstat bands
cl <- colorRampPalette(c('darkorchid3','light blue','lightpink4'))(100) 
plot(dvi2011)

#EXERCISE: plot DVI for 1988
dvi1988 <- p224r63_1988_masked$B4_sre - p224r63_1988_masked$B3_sre
cl <- colorRampPalette(c('darkorchid3','light blue','lightpink4'))(100) 
plot(dvi1988)

#Let's see the difference between the two years
diff <- dvi2011 - dvi1988

#See the effects of changing the scale, resolution
#resolution = dimension of pixels = resempling
p224r63_2011res <- aggregate( p224r63_2011_masked, fact=10) # "fact" increase the pixel size * 10
p224r63_2011res100 <- aggregate( p224r63_2011_masked, fact=100) #increase pixel size * 100

#plot the regular and enhanced images all together in column 
par(mfrow=c(3,1))
plotRGB(p224r63_2011_masked, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res, r=4, g=3, b=2, stretch="Lin")
plotRGB(p224r63_2011res100, r=4, g=3, b=2, stretch="Lin")

############################################################
############################################################
############################################################

########## 7. R code for ecosystem function
# R code to see biomass over the world and calculate changes in ecosystem functions (energy, chemical cycling and proxies)

install.packages("rasterdiv") #this package provides functions to calculate indices of diversity on raster data
install.packages("rasterVis") #package for visualization
library(rasterdiv) 
library(rasterVis) 

data(copNDVI) #Copernicus Long-term database
plot(copNDVI)

copNDVI <- reclassify(copNDVI, cbind(253:255, NA)) # removing water-based colours -> reclassify: function that reclassify groups of values to other values
levelplot(copNDVI) #plot the image as a level plot               

copNDVI10 <- aggregate (copNDVI, fact=10) # aggregating 10 pixel in 1 
levelplot(copNDVI10) 
copNDVI100 <- aggregate (copNDVI, fact=100) #exageratig the aggregation
levelplot(copNDVI100) 

###### Let's continue
setwd("C:/LAB/")
library(raster)

defor1 <- brick("defor1_.jpg")
defor2 <- brick("defor2_.jpg")

#REMEMBER:
# Band 1: NI
# Band 2: red
# Band 3: green 

#plot the images together
par(mfrow=c(1,2))
plotRGB(defor1, r=1, g=2, b=3, stretch="Lin")
plotRGB(defor2, r=1, g=2, b=3, stretch="Lin")

#Let's calculate DVI for defor1 and defor2
dvi1 <- defor1$defor1_.1 - defor1$defor1_.2 # we use $ because we haven't attached defor1 to R
dvi2 <- defor2$defor2_.1 - defor2$defor2_.2

#Plot the two graphs for DVI with a specific colour scheme
cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100) # specifying the color scheme
par(mfrow=c(1,2))
plot(dvi1, col=cl)
plot(dvi2, col=cl)

dev.off() 

difdvi <- dvi1 - dvi2
cld <- colorRampPalette(c('blue','white','red'))(100) 
plot(difdvi, col=cld)

hist(difdvi) #create an histogram 

############################################################
############################################################
############################################################

########## 8. R code for PCA remore sensing

setwd("C:/LAB")
library(raster)
library(RStoolbox)
library(ggplot2)

p224r63_2011 <- brick("p224r63_2011_masked.grd") # import the images 

#B1: blue
#B2: green
#B3: red
#B4: NIR
#B5: SWIR (Short-Waved InfraRed)
#B6: thermal infrared 
#B7: SWIR
#B8: panchromatic

#Plot the image with RGB
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin") 
ggRGB(p224r63_2011,5,4,3) # "ggRGB" calculates RGB color composite raster for plotting with ggplot2

#Now plot the 1988 image
p224r63_1988 <- brick("p224r63_1988_masked.grd")
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
ggRGB(p224r63_1988,5,4,3)

#Plot both images together (2011 and 1988)
par(mfrow=c(1,2))
plotRGB(p224r63_1988, r=5, g=4, b=3, stretch="Lin")
plotRGB(p224r63_2011, r=5, g=4, b=3, stretch="Lin")

names(p224r63_2011) #see the bands
plot(p224r63_2011$B1_sre,p224r63_2011$B1_sre) # to see if these bands are correlated 

###### PCA -> summarize the information in large data tables by means of a smaller set of “summary indices” that can be more easily visualized and analyzed
#decrease the resolution 
p224r63_2011_res <- aggregate(p224r63_2011, fact=10)
p224r63_2011_pca <- rasterPCA(p224r63_2011_res) #Calculates PCA for RasterBricks or RasterStacks and returns a RasterBrick with multiple layers of PCA scores.

plot(p224r63_2011_pca$map) #in this way we have the default color palette

#Let's create a gray plot
cl <- colorRampPAlette(c("dark grey","grey","light grey"))(100)
plot(p224r63_2011_pca$map, color=cl) 

#See the results of the function we have applied
summary(p224r63_2011_pca$model)

#Let's plot the first 3 components of the 2011 image
pairs(p224r63_2011)
plotRGB(p224r63_2011_pca$map, r=1, g=2, b=3, stretch="Lin")

###### Now let's do the same of the 1988 image
p224r63_1988_res <- aggregate(p224r63_1988, fact=10)
p224r63_1988_pca <- rasterPCA(p224r63_1988_res) 
plot(p224r63_1988_pca$map, color=cl) 
summary(p224r63_1988_pca$model)
pairs(p224r63_1988)

#Calculate the difference in PCA between the two images
diffpca <- p224r63_2011_pca$map - p224r63_1988_pca$map
plot(diffpca)
cldif <-colorRampPalette(c('blue','black','yellow'))(100) 
plot(diffpca$PC1, color=cldif)

############################################################
############################################################
############################################################

########## 9. R code for radiance

library(raster)
toy <- raster(ncl=2, nrow=2, xmn=1, xmx=2, ymn=1, ymx=2)
values(toy) <- c(1.13, 1.44,1.55,3.4)
plot(toy)
text(toy,digits=2)

# 2 bits --> 2^2 = 4 values 
toy2bit<- stretch(toy,minv=0, maxv=3) #stretch changes the range of values -> minimum of 0 to a maximum of 3

storage.mode(toy2bits[]) = "integer"  #to ensure that we are going to use integer value 
plot(toy2bits)
text(toy2bits, digits=2)    

# 4 bits --> 2^4 =16 values 
toy4bit<- stretch(toy,minv=0, maxv=15)
storage.mode(toy4bits[]) = "integer"
plot(toy4bits)
text(toy4bits, digits=2) 

# 8 bits --> 2^8 =256 values 
toy8bit<- stretch(toy,minv=0, maxv=255)
storage.mode(toy8bits[]) = "integer"
plot(toy8bits)
text(toy8bits, digits=2)

#plot all together
par(mfrow=c(1,4))
  plot(toy)
text(toy,digits=2)
  plot(toy2bits)
text(toy2bits, digits=2)
  plot(toy4bits)
text(toy4bits, digits=2) 
  plot(toy8bits)
text(toy8bits, digits=2)

############################################################
############################################################
############################################################

########## 10. R code for faPAR






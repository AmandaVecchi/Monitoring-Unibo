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

setwd("C:/LAB/")
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
#faPAR = fraction of the solar radiation absorbed by living leaves 

setwd("C:/LAB/")
library(raster)
library(rasterVis)
library(rasterdiv)
library(sf)

plot(copNDVI) # Copernicus NDVI from the rasterdiv package

#Let's reclassify the data and remove water from the analysis
copNDVI <- reclassify(copNDVI, cbind(253:255, NA)) # "cbind" removes data, what value to substitute to the removed ones
# NA = no data
levelplot(copNDVI)

faPAR10 <- raster("farPAR10.tif") # import the image
levelplot(farPAR10)

#Save the plot as .PDF
pfd("copNDVI.pdf")
levelplot(coNDVI)
dev.off()

pdf("faPAR.pdf")
levelplot(faPAR10)
dev.off()

#######Let's continue
setwd("C:/LAB/") 
load("faPAR.RData")
library(raster)
library(rasterdiv)
library (rasterVis)

#the original faPAR from Copernicus is 2GB. Let's see how much space is needed for an 8-bit set
writeRaster(copNDVI, "copNDVI.tif") # 5.3 MB

levelplot(faPAR10) 

########## Regression model between faPAR and NDVI

erosion <- c(12,14,16,24,26,40, 55,67) #ex. amount of erosion in a certain area
hm <- c(30,100,150,200,260,340,,460,600) #ex. amount of heavy metals

plot(erosion, hm, col="red", pch=19, xlab="Erosion", ylab="Heavy metals") #xlab and ylab set the lables for x and y axis

#Create a linear model!!!
model1 <- lm(hm ~ erosion)
summary(model1)
abline(model1) #function that creates regression lines on plots
#We create a regression line that relates erosion and number of heavy metals

## faPAR vs NDVI model
setwd("C:/LAB/")
library(raster)
faPAR10 <- raster("farPAR10.tif")
plot(faPAR10)

plot(copNDVI)
copNDVI <- reclassify(copNDVI, cbind(253:255, NA), right=TRUE)

#We want to see how the two variables are related
# RANDOM SAMPLES
random.points <- function(x,n)  # x is the raster file, n is the number of the random points
{
lin <- rasterToContour(is.na(x))
pol <- as(st_union(st_polygonize(st_as_sf(lin))), 'Spatial') 
pts <- spsample(pol[1,], n, type = 'random')
}

pts <-  random.points(faPAR10, 1000) #ex. we select 1000 points from faPAR10
plot(faPAR10)
points(pts,col="red",pch=19)
 
#Let's extract points from a raster
copNDVIp <-extract (copNDVI, pts)
faPAR10p <-extract (faPAR10, pts)

#Create a model for PHOTOSYNTHESIS vs BIOMASS
model2 <- lm(faPAR10p ~ copNDVIp)
plot(copNDVIp, faPAR10p, col="green", xlab="Biomass", ylab="Photosynthesis")
abline(model2, col="red")

############################################################
############################################################
############################################################

########## 11. R code EVB (measure the standard deviation from a satellite image)

setwd("C:/LAB/")

library(raster)
library(Rstoolbox) #for PCA

snt <- brick("snt_r10.tif")
plot(snt)

#REMEMBER:
#B1 = blue
#B2 = green
#B3 = red
#B4 = NIR
#RGB -> R = 3, G = 2, B = 1

plotRGB(snt, 3,2,1, stretch="lin") # plot the image with RGB system, visible colours 
plotRGB(snt, 4,3,2, stretch="lin") # plotting NIR in top of red = vegetation coloured in red

pairs(snt) #Scatterplot matrices

###### PCA analysis
sntpca <- rasterPCA(snt) # "rasterPCA" calculates R-mode PCA for raster images
sntpca #see the result
summary(sntpca$model)
plot(sntpca$map) 

# Make a RGB plot
plotRGB(sntpca$map, 1,2,3, stretch="lin") 

#Let's calculate che standard deviation and create a moving window
window <- matrix(1, nrow=5, ncol=5) #"matrix" function is used to create a matrix defined by the arguments nrow and ncol
window #here we have the moving window

#"focal" function calculates values for the neighborhood of focal cells component
sd_snt <- focal(sntpca$map$PC1, w=window, fun=sd) 
# "w" states the moving window we want to use
# "fun" states the function to be applied. sd = standard deviation
cl <- colorRampPalette(c("dark blue", "green", "orange", "red"))(100)
plot(sd_snt, col=cl)

par(mfrow=c(1,2))
plotRGB(snt,4,3,2, stretch="lin", main="original image") 
plot(sd_snt, col=cl, main="diversity")

######Let's continue
#Focus on Cladonia stellaris species

library(raster)
library(RStoolbox)

setwd("C:/LAB/")

clad <- brick("cladonia_stellaris_calaita.JPG")  #Importing the image related to C. stellaris

#Create a moving window of a matrix of 3 by 3 pixels (number 1 doesn't impact the calculation) 
window <- matrix(1, nrow = 3, ncol = 3)
window

#PCA analysis for C. stellaris
cladpca <- rasterPCA(clad)
cladpca #see the output of the applied function
summary(cladpca$model) #variance is 0.98 = 98% of the first component 
plotRGB(cladpca$map, 1,2,3, stretch="lin") #plotting the results

#Set the moving window again
window <- matrix(1, nrow=5, ncol=5)
window

# "focal" calculates values for neighborhood of focal cells
sd_clad <- focal(cladpca$map$PC1, w=window, fun=sd) 

#Apply the aggregate function
PC1_agg <- aggregate(cladpca$map$PC1, fact=10)
sd_clad_agg <- focal(PC1_agg, w=window, fun=sd) 

#Let's see both graphs together to see the differencies 
par(mfrow=c(1,2))  
cl <- colorRampPalette(c('yellow','violet','black'))(100) 
plot(sd_clad,col=cl) 
plot(sd_clad_agg,col=cl)  #cladonia set aggregated 

#Let's plot the calculation 
par(mfrow=c(1,2)) 
cl <- colorRampPalette(c('yellow','violet','black'))(100) 
plotRGB(clad, 1,2,3, stretch="lin")
plot(sd_clad, col=cl)

dev.off()

############################################################
############################################################
############################################################

########## 12. R code NO2

setwd("C:/LAB/NO2/")
install.packages("ncdf4") 
library(ncdf4) #ncdf4 is a package to read all the netCDF files. All Copernicus data use this extension

#EXERCISE:import all of the NO2 data in R using the lapply function 
#We want to import a set of files, so we create a list spcifying all files we want to import
rlist <- list.files(pattern="EN")  # "pattern" = states on which base R has to select the imported data: they all contain in the name "EN" 

import <- lapply(rlist, raster) #importing all the file selected with rlist 

EN <- stack(import) #"stack" function is used to create a multitemporal image 
cl <-colorRampPalette(c('red','orange','yellow'))(100) 
plot(EN,col=cl) #plot the images all together 

#Let's see the difference from the first and last day 
par(mfrow=c(1,2))
plot(EN$EN_0001, col=cl)
plot(EN$EN_0013, col=cl)

#Plot with RGB the first, seventh and last image
plotRGB(EN, r=1, g=7, b=13,stretch="lin") 

#Let's see the difference between the final and initial day
dif <- EN$EN_0013 - EN$EN_0001
cld <-colorRampPalette(c('blue','white','red'))(100) #red means high differences, blue means lower difference
plot(dif,col=cld) 

#Quantitative measure of the decreasing of NO2
 #Let's create a boxplot to see the relative differencies
boxplot(EN)
boxplot(EN,outline=F) #vertical boxplot with the outlines removed 
boxplot(EN,outline=F, horizontal=T) #horizontal boxplot
boxplot(EN,outline=F, horizontal=T, axes=T)

plot(EN$EN_0001, EN$EN_0013)

# Comparing each pixels of 2 stituations plotting the 2 images in order to see if NO2 decreases o increases in each pixel
plot(EN$EN_0001, EN$EN_0013) 
abline(0,1,col="red") #most of the point are under the line, that means that the NO2 decreases

############################################################
############################################################
############################################################

########## 13. R code snow

setwd("C:/LAB/SNOW/")
library(ncdf4)
library(raster)

#Let's import the required image
snowmay <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc")

cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 

#EXERCISE: plot snow cover with the cl palette
plot(snowmay,col=cl) 

#Now let's import the snow related data
snow2000r <- raster("snow2000r.tif")
snow2005r <- raster("snow2005r.tif")
snow2010r <- raster("snow2010r.tif")
snow2015r <- raster("snow2015r.tif")
snow2020r <- raster("snow2020r.tif")

#Let's look at all the data toether represented with the same colour palette
par(mfrow=c(2,3)) #multiframe row 
plot(snow2000r, col=cl)
plot(snow2005r, col=cl)
plot(snow2010r, col=cl)
plot(snow2015r, col=cl)
plot(snow2020r, col=cl)

####FASTER version to mport all the required images
# rlist <- list.files(pattern="snow")  #they all contain the word "snow"--> same pattern 
# import<- lapply(rlist, raster) 
# snow.multitemp <- stack(import) # "stack" function allows us to create a multitemporal image 
plot(snow.multitemp,col=cl) # we have all the plot together without using the par function 

###Let's make a prediction
source("prediction.r") #we use a script in R 

# Since the code needs time, we can download predicted.snow.2025.norm.tif from iol in the Data and upload it into R
predicted.snow.2025.norm <- raster("predicted.snow.2025.norm.tif")
plot(predicted.snow.2025.norm, col=cl)

######Let's continue

setwd("C:/LAB/SNOW")

#EXERCISE: import the snow cover images all together
rlist <- list.files(pattern="snow")
import <- lapply(rlist, raster)
snow.multitemp <- stack(import)
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100) 
plot(snow.multitemp, col=cl)

#Load the prediction file (we did it the previous time already) 
prediction <- raster("predicted.2025.norm.tif")
plot(prediction, col=cl) #snow cover will be present only in the northern part of the world 

#Now export the output we obtained
writeRaster(prediction, "final.tif") 
#in the folder "SNOW" we have now have a file named "final.tif

final.stack <- stack(snow.multitemp, prediction)
plot(final.stack, col=cl) #plot all the images together 

#export the R graph as a pdf
pdf("my_final_graph.pdf")
plot(final.stack, col=cl)
dev.off()

#export the R graph as a png
png(("my_final_graph.pdf")
plot(final.stack, col=cl)
dev.off()

############################################################
############################################################
############################################################

########## 14. R code "CROP"

setwd("C:/LAB/") 
    
library(raster)
library(ncdf4)
 
snow <- raster("c_gls_SCE_202005260000_NHEMI_VIIRS_V1.0.1.nc")
cl <- colorRampPalette(c('darkblue','blue','light blue'))(100)
plot(snow, col=cl)

#Let's have a look only to Italy
ext <- c(0, 20, 35, 50)

#Zoom on the desired area
zoom(snow, ext=ext)

#Crop the image with the focus on the area that we want 
crop(snow, ext)
snowitaly <- crop(snow, ext)

#to obtain a rectangular image
zoom(snow, ext=drawExtent())

############################################################
############################################################
############################################################

########## 15. R code interpolation

setwd("C:/LAB/")

library(spatstat)
inp <- read.table("dati_plot55_LAST3.csv", sep=";", head=T)
head(inp)  

attach(inp) #attach our dataframe (inp)
#plot(inp$X, inp$Y) if we do not use the "attach" function
plot(X,Y)
summary(inp)

#Let's make a planar point pattern
inppp <- ppp(x=X, y=Y, c(716000,718000),c(4859000,4861000)) #Stating the 2 variables and the range "c()" for them

names(inp) #Let's have a look at the names of the variables

marks(inppp) <- Canopy.cov #Extracting the marks

canopy <- Smooth(inppp) #we are going to interpolate the data 
plot(canopy) #plot smooth
points(inppp,col="green")

##Lichens cover
marks(inppp)<- cop.lich.mean  #lichens on the trees and absence of lichens reflects high pollution 
lichs<- Smooth(inppp)
plot(lichs)
points(inppp) #We see that the higher amount of lichens is in the north, less polluted area

#Plot the images together
par(mfrow=c(1,2))
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)
plot(Canopy.cov, cop.lich.mean, col="red", pch=19, cex=2)

### Psammophilous vegetation
inp.psam <- read.table("dati_psammofile.csv", sep=";", head=T)
attach(inp.psam) 
head(inp.psam)

plot(E,N) 
inp.psam.ppp <- ppp(x=E,y=N,c(356450,372240),c(5059800,5064150))

marks(inp.psam.ppp)<- C_org # measure carbon in the soil 
C <- Smooth(inp.psam.ppp) #C = carbon
plot(C) 
points(inp.psam.ppp)

############################################################
############################################################
############################################################

########## 16. R code SDM = Species Distribution Modelling
    
install.packages("sdm") 
install.packages("rgdal")
library(sdm)
library(raster) #used for geological variables: predictors  (shp)
library(rgdal) #used for species

#Let's import species data
file<-system.file("external/species.shp", package="sdm") 
species<-shapefile(file)
species #too see all the information about the species 
    
species$Occurrence #occurrence is 0 (= absent) or 1 ( = present)
plot(species[species$Occurrence == 1,],col="blue",pch=16)
points(species[species$Occurrence == 0,],col="red",pch=16)
    
#Environmental variables
path<-system.file("external", package="sdm") 
lst<-list.files(path=path,pattern="asc$",full.names = T) 
lst

#Plotting all the variables inside the stack
preds<-stack(lst)
cl<-colorRampPalette(c("blue","orange","red","yellow")) (100)
plot(preds, col=cl)

plot(preds$elevation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$temperature, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$precipitation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$vegetation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

#Let's make the model
d<-sdmData(train=species, predictors=preds) #training set=insitu data(species set make previously), predictors=where the sp are predicted to be (stack of all of the variables) 
d #to see what there is inside 
m1<-sdm(Occurrence~elevation+precipitation+temperature+vegetation, data=d, methods="glm") 
p1<-predict(m1, newdata=preds)
plot(p1, col=cl)
points(species[species$Occurrence == 1,], pch=16)
s1 <- stack(preds,p1)
plot(s1, col=cl)
    
############################################################    
############################################################    
############################################################    
    
########## 17. R code exam project
    
install.packages("raster")
install.packages("rasterVis")
install.packages("rasterdiv")
install.packages("rgdal")
install.packages("ggplot2")
install.packages("RStoolbox")
install.packages("gridExtra")
    
library(raster)
library(rasterVis)
library(rasterdiv)
library(rgdal)
library(ggplot2)
library(RStoolbox)
library(gridExtra)
    
#Upload image from 2016
setwd("C:/empedrado/prefire/")
rlist2016 <- list.files( pattern = "2016")
rlist2016

#B2 = BLUE
#B3 = GREEN
#B4 = RED
#B5 = NIR
#B6 = SWIR1
#B7 = SWIR2

import_image2016 <- lapply(rlist2016, raster)
image2016 <- stack(import_image2016)
image2016

#crop the image
ext <- c(-8100000, -8000000, -4325000, -4275000)
zoom(image2016, ext=ext)
image2016_crop <- crop(image2016, ext)
plot(image2016_crop)

#plot the image with RGB colours 
RGB2016 <- ggRGB(image2016_crop, r=3, g=2, b=1, stretch ="Lin")
plot(RGB2016)
    
visible2016 <- print(RGB2016 + ggtitle("Visible - Empedrado, 2016")+ theme (axis.title.x = element_blank(),axis.title.y = element_blank()))    
    
setwd("C:/empedrado/postfire/")
rlist2019 <- list.files(pattern = "2019")
rlist2019
    
import_image2019 <- lapply(rlist2019, raster)
image2019 <- stack(import_image2019)
image2019

#crop the image
ext <- c(-8100000, -8000000, -4325000, -4275000)
zoom(image2019, ext=ext)
image2019_crop <- crop(image2019, ext)
plot(image2019_crop)

RGB2019 <- ggRGB(image2019_crop, r=4, g=3, b=2, stretch ="Lin")
plot(RGB2019) 
    
visible2019 <- print(RGB2019 + ggtitle("Visible - Empedrado, 2019")+ theme (axis.title.x = element_blank(),axis.title.y = element_blank()))    

grid.arrange(visible2016, visible2019, nrow =1)

#plot both RGB images together
par(mfrow = c(1,2) , oma=c(0,0,2,0))
plotRGB(image2016_crop, r=4, g=3, b=2, stretch ="lin")
plotRGB(image2019_crop, r=4, g=3, b=2, stretch ="lin")
mtext("Visible Empedrado 2016 vs Empedrado 2019", outer = TRUE, cex = 1.5)

par(mfrow = c(1,2), oma = c(0,0,2,0))
plotRGB(image2016_crop,  r=4, g=3, b=2, stretch ="hist")
plotRGB(image2019_crop, r=4, g=3, b=2, stretch ="hist")
mtext("Visible enhanced Empedrado 2016 vs Empedrado 2019", outer=TRUE, cex =1.5)
    
#Assess Plant density and health
false2016 <- ggRGB(image2016_crop, r=5, g=4, b=3,stretch = "Lin") 
false2019 <- ggRGB(image2019_crop, r=5, g=4, b=3,stretch = "Lin") 

par(mfrow = c(1,2) , oma=c(0,0,2,0))
plot(false2016, main = "Fase colour - 2016")
plot(false2019, main = "False colour - 2019")
mtext("False colour: Empedrado 2016 vs Empedrado 2019", outer = TRUE, cex = 1.5)
    
#CALCULATE NDVI
#NDVI = (NIR – RED) / (NIR+RED)
NDVI2016 <- (image2016_crop[[5]] - image2016_crop[[4]]) / (image2016_crop[[5]] + image2016_crop[[4]])
NDVI2016
plot(NDVI2016, main = "NDVI - Empedrado, 2016")

NDVI2019 <- (image2019_crop[[5]] - image2019_crop[[4]]) / (image2019_crop[[5]] + image2019_crop[[4]])
NDVI2019
plot(NDVI2019, main = "NDVI - Empedrado, 2019")

#PLOT BOTH NDVI IMAGES TOGETHER
par(mfrow = c(1,2), oma=c(0,0,2,0))
plot(NDVI2016, main = "NDVI 2016") 
plot(NDVI2019, main = "NDVI 2019") 
mtext("NDVI comparison: 2016 vs 2019", outer=TRUE, cex =1.5)

#Difference between the two years
diffNDVI <- (NDVI2016 - NDVI2019)
plot(diffNDVI)
plot(diffNDVI, main = "NDVI difference between 2016 and 2019")

cld <- colorRampPalette(c('blue','white','red'))(100)
plot(diffNDVI, col = cld, main = "NDVI difference between 2016 and 2019")

# remove cells where NDVI < 0,4
NDVI2016_mod <- reclassify(NDVI2016, cbind(-Inf, 0.4, NA))
NDVI2019_mod <- reclassify(NDVI2019, cbind(-Inf, 0.4, NA))
     
par(mfrow=c(1,2))
plot(NDVI2016_mod, main="Vegetation 2016", axes=FALSE)
plot(NDVI2019_mod, main="Vegetation 2019", axes=FALSE)

#VISUALIZATION USING HISTOGRAMS
hist_NDVI2016 <- hist(NDVI2016, main = "Distribution of NDVI values - 2016", xlab = "NDVI", ylab = "Frequency", breaks = 50)
hist_NDVI2016 <- hist(NDVI2016, main = "Distribution of NDVI values - 2016", xlab = "NDVI", ylab = "Frequency", breaks = 50, col = "light pink",border = "black")
     
hist_NDVI2019 <- hist(NDVI2019, main = "Distribution of NDVI values - 2019", xlab = "NDVI", ylab = "Frequency", breaks = 50)
hist_NDVI2019 <- hist(NDVI2019, main = "Distribution of NDVI values - 2019", xlab = "NDVI", ylab = "Frequency", breaks = 50, col = "light blue",border = "black")
     
par(mfrow = c(1,2), oma=c(0,0,2,0))
plot(hist_NDVI2016, main = "NDVI 2016", col = "light pink",border = "black") 
plot(hist_NDVI2019, main = "NDVI 2019", col = "light blue",border = "black") 
mtext("Comparison of distribution of NDVI values: 2016 vs 2019", outer=TRUE, cex =1.5)
     
col2rgb("lightblue")
    # R 173
    # G 216
    # B 230
col2rgb("pink")
    # R 255
    # G 192
    # B 203
col2rgb(c("lightblue", "lightgreen", "pink"))

mycol <- rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50") #make the color transparent by 50%
c1 <- rgb(173,216,230, max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "pink") 

par(mai=rep(0.5, 4)) #set the margins for the image
layout(matrix(c(1,1,2,2,0,3,3,0), ncol = 4, byrow = TRUE)) #devide the plotting space
plot(hist_NDVI2016, col=c2, main="NDVI 2016", xlab = "NDVI")
plot(hist_NDVI2019, col=c1, main="NDVI 2019", xlab = "NDVI")
plot(hist_NDVI2016, col = c2, xlim = c(-0.5, 1), main="Comparison between NDVI",xlab = "NDVI" )
plot(hist_NDVI2019, add = TRUE, col = c1)

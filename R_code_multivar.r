# R code for multivariate analysis

#Setting the working diretory and istalling + stating we are using the "vegan" package
#The vegan package provides tools for descriptive community ecology; it contains tools for diversity analysis on vegetation
setwd("C:/LAB/")
install.packages("vegan")
library(vegan) 

#Download dataset from IOL and import them in R
#We are giving a name to the dataset, assigning a name to the read.table function (<-)
biomes <-read.table("biomes.csv", header=T, sep=",") #sep is an argument that declares which is the separator of the lines in the dataset
head(biomes) #view(biomes) 

# Multivariate analysis
#The "decorana" function performs detrended correspondence analysis (=see data in 2D)
#DCA is a technique used to find the main factors or gradients in large, species-rich but usually sparse data matrices 
multivar <- decorana(biomes) #We are also giving a name to the correspondence analysis we are performing
plot(multivar)

#See how much data we are loosing when transforming the data in a 2D vision
multivar

plot(multivar)
biomes_types <-read.table("biomes_types.csv", header=T, sep=",") #We are making use of a new dataset
head(biomes_types)
attach(biomes_types) #attach our dataset

#The "ordiellipse" function draws ellipses connecting all point belonging to one biome
ordiellipse(multivar, type, col=1:4, kind="ehull", lwd=3) 
#We declare the dataset and the columns we want to base our ellipses on
#col= colours to distinguish different biomes
#kind = "ehull" = draws ellipsoid that enclose all points in the group (ehull)
#lwd = dimension of line of ellipse
#Or also col=c("green", "blue", "black", "red")

ordispider(multivar, type, col=1:4, label=TRUE) #In this case we are creating a spider diagram 








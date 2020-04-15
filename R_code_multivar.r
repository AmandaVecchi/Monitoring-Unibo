# R code for multivariate analysis

library(vegan) #package already installed

#download dataset from IOL and import them in R
#give a name to dataset + assign name to a function (<-) + read.table to import from outside (name) + header = true bc first line is title + separator is the comma, so we state it
biomes <-read.table("biomes.csv", header=T, sep=",")

#look at the dataset
head(biomes) #view(biomes) # biomes

# Multivariate analysis
# detrended correspondence analysis (=see data in 2D) --> decorana
multivar <- decorana(biomes)
plot(multivar)

#see how much data we are loosing when transforming in a 2D vision
multivar # eigenvaues= perception of the whole system in all 4 dimensions(expressed in percentage)

plot(multivar)
biomes_types <-read.table("biomes_types.csv", header=T, sep=",")
head(biomes_types)

#link the points to see if different biomes can be seen in ou graph
attach(biomes_types) #attach our dataset
#draw ellipse connecting all point of one biome
ordiellipse(multivar, type, col=1:4, kind="ehull", lwd=3) #declare the colums to make group + col= colours to distinguish different biomes + kind of raggruppation +  dimension of line of ellipse
# col=c("green", "blue", "black", "red")

ordispider(multivar, type, col=1:4, label=TRUE) #species not included in the ellipse are those that cannot be included in the ellipse bc we are looking at our dataset in only 2D out of 4D










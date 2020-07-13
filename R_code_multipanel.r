## Multipanel in R: second lecture of Monitoring Ecosystems

install.packages("sp")
install.packages("GGally)
#GGally is an extension of ggplot2; it adds functions to reduce the complexity of combining geometric objects with transformed data

#library = tells R we want to use the package installed
library(sp)
#require(sp) will also do the job
library(GGally)

data(meuse) 
attach(meuse)

#names = function to get the name of all variables in the dataset
names(meuse) 
head(meuse) #head = function that shows only the first six lines of the dataset

#Plot cadmium vs zinc
plot(cadmium,zinc)
plot(cadmium, zinc, pch=15, col="red", cex=2) 
#pch = argument used to specify point shapes
#col = declare the color of the point plotted
#cex = number indicating the amount by which plotting symbols should be scaled relative to the default. 2 = all points enhanced by 100%

#pairs = creates matrix of all scatterplots of the dataset
pairs(meuse) #creates all the possible pairwise plots of the dataset

#error = "size of the margins is too large". Solution = reshape with mouse the graph window and relaunch the code

#Create pairs scatterplots considering only some variables:cadmium, copper, lead, zinc
pairs(~ cadmium + copper + lead + zinc, data=meuse)
pairs(meuse[,3:6]) #start from column number 3 to column number 6 

#Prettify this graph
pairs(meuse[,3:6], pch=15,col="blue") #Change the symbol and colour
pairs(meuse[,3:6], pch=15,col="blue", cex=0.5) #Change the dimension of the symbol; 505 smaller

#Add info to the scatter plot with GGally (install a new package, but at beginning of code)
#use the library (at beginning of code) 
ggpairs(meuse[,3:6])




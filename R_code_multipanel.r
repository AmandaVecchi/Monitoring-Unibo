## Multipanel in R

install.packages("sp")
install.packages("GGally")
#GGally is an extension of ggplot2; it adds functions to reduce the complexity of combining geometric objects with transformed data

library(sp) #require(sp) 
library(GGally)

data(meuse) 
attach(meuse)

#The "names" function allows usto get the name of all variables in the dataset
names(meuse) 
head(meuse) 

#Plot cadmium vs zinc
plot(cadmium,zinc)
plot(cadmium, zinc, pch=15, col="red", cex=2) 
#pch = argument used to specify point shapes
#col = argument used to declare the color of the point plotted
#cex = number indicating the amount by which plotting symbols should be scaled relative to the default

#The "pairs" function creates a matrix of all scatterplots of the dataset
pairs(meuse) #creates all the possible pairwise plots of the dataset

#IF we get error = "size of the margins is too large", the solution is to reshape with our mouse the graph window and relaunch the "pairs" code

#Create pairs scatterplots considering only some variables:cadmium, copper, lead, zinc
#The tilde operator states which are the dependent variables
pairs(~ cadmium + copper + lead + zinc, data=meuse) #in this case we also state the dataset in which the variables are contained
pairs(meuse[,3:6]) #using the square brackets tell R to start take all variables starting from column number 3 to column number 6 in the dataset we are using

#Let's prettify this graph
pairs(meuse[,3:6], pch=15,col="blue") 
pairs(meuse[,3:6], pch=15,col="blue", cex=0.5) #cex = 0.5 diminishes point size of 50%

#Now let's add info to the scatter plot with the GGally package (install a new package, but at beginning of code)
#The ggpairs() function of the GGally package allows to build a scatterplot matrix
ggpairs(meuse[,3:6])




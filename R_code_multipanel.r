## Multipanel in R: second lecture of Monitoring Ecosystems

install.packages("sp")
install.packages("GGally)

#tell R we want to use the package installed
library(sp) # require(sp) will also do the job
library(GGally)

#choose the dataset we want to use
data(meuse) # there is an available dataset named meuse

#use the dataset 
attach(meuse)

#see the names of the variables in the dataset
names(meuse) #see all names 
head(meuse) #shows only the first six lines 

#plot cadmium vs zinc
plot(cadmium,zinc)
plot(cadmium,zinc,pch=15,col="red",cex=2) # change symbol, how big it is and colour

#make all the possible pairwise plots of the dataset
pairs(meuse) 

#error = "size of the margins is too large" --> reshae with mouse the graph window and relaunch the code

#consider a lower amount of variables --> cadmium, copper, lead, zinc
pairs(~ cadmium + copper + lead + zinc, data=meuse)
pairs(meuse[,3:6]) #start from column number 3 to column number 6 

#prettify this graph
pairs(meuse[,3:6], pch=15,col="blue") #changed the symbol and colour
pairs(meuse[,3:6], pch=15,col="blue", cex=0.5) #changed also the dimension of the symbol

#add info to the scatter plot --> istall a new package, but at beginning of code
#use the library --> at beginning of code to download the package
#GGally will prettify the plot
ggpairs(meuse[,3:6])




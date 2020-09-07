## Our first R code

#The "install.packages" function downloads and install packages from CRAN repositories or from local files.
install.packages("sp") #sp is a package that allows us to manipulate and work on spatial datasets
#The "library" or "require" function specifies to R we want to use the package we have installed
library(sp)

#The "data" function allows us to upload specified data sets
data(meuse) #meuse is a dataset containing locations and topsoil heavy metal concentrations

#Let's see how the meuse dataset is structured
meuse

#The "head" function allows us to see the first six rows of the data set we are using
head(meuse)

#Let's plot two variables together to see if they are related
#See if zinc concentration is related to that of copper 

#The "attach" function connects the database to the R search path. With this function the objects (in this case heavy metal concentrations) in the database can be accessed by giving their names
attach(meuse)

#The "plot" function plots on a generic x;y graph two different variables
plot(zinc, copper)
plot(zinc, copper,col="green") #col = argument to state in which colour to plot the points in the graph
plot(zinc, copper,col="green",pch=19) #pch = argument used to specify point shapes
plot(zinc, copper,col="green",pch=19,cex=2) #cex = argument used to specify character exaggeration 

#cex = 0.5 diminishes the size of the point of 50%
#cex= 2 enhances the size of the point of 200%
#cex = 0.7 decreases the seize of the poin of 70%
#in general to mininish the size of the pont we put 0 < cex < 1

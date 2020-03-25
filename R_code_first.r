install.packages("sp")

library(sp)
data(meuse)

#Let's see how the meuse dataset is structured:
meuse

#Let's look at the first row of the set. Only see the first part of the dataset, six lines
head(meuse)

#Let's plot two variables together to see if they are related
#See if zinc concentration is related to that of copper --> attach a dataset to R first, then plot the two variables
attach(meuse)
plot(zinc, copper)
plot(zinc, copper,col="green")
#change symbol of plot --> pch (=symbol) and the number of the symbol
plot(zinc, copper,col="green",pch=19)
#make bigger point --> cex = character exhageration --> 2= double the symbol. To decrease put number between 0 and 1 
plot(zinc, copper,col="green",pch=19,cex=2)

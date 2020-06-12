#species distribution modelling 

install.packages("sdm")
library(sdm)
library(raster) #ecological variables = predictors = ecological variables used to predict species distribution
# install.packages("rgdal")
library(rgdal) #species, 

#import species data
file <- system.file("external/species.shp", package="sdm")
species <- shapefile(file) #shapefile is the extent of files downloaded with the installing of the sdm package
plot(species) #points represent if the specis is present 
species #see the info of the dataset

plot(species[species$Occurrence == 1,],col='blue',pch=16)
# [] condition --> species occurrence has to be 1; put == to look only at data with a value equals to something
 
 points((species[species$Occurrence == 0,],col='red',pch=16) 
 #with point function we add point of non-occurrence to previous map. If we use "plot function we are deleting the previous map
 
 #import the predictors
 path <- system.file("external", packages="sdm") #path towards the external folder
 
 #create a list of all predictors
 #pattern = select files that have in the name this part
 lst <- list.files(path=path,pattern='asc$',full.names = T) 
 lst
 
 #make a stack of all data
preds <- stack(lst)

cl <- colorRampPalette(c('blue','orange','red','yellow')) (100)
plot(preds, col=cl)

plot(preds$elevation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$temperature, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$precipitation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

plot(preds$vegetation, col=cl)
points(species[species$Occurrence == 1,], pch=16)

#put all info together and make a model
#explain the data to make the model
d <- sdmData(train=species, predictors=preds)
d

m1 <- sdm(Occurrence ~ elevation + precipitation + temperature + vegetation, data=d, methods = "glm")
p1 <- predict(m1, newdata=preds)

plot(p1, col=cl)
points(species[species$Occurrence == 1,], pch=16)

s1 <- stack(preds, p1)
plot(s1, col=cl)
































 





































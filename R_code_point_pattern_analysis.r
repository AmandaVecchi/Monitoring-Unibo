## Point pattern analysis: density maps

install.packages("spatstat")
library(spatstat)

#attach our dataset, the one we are going to use
attach(covid)

#what are the coordinates of the datases and what are the extensions of these cooridnates
head(covid)
covids <- ppp(lon, lat, c(-180,180), c(-90,90))  #ppp= planar point pattern, specify the variables and the range

#density map
d <- density(covids)

#show density map
plot(d)

#put points on the density map
points(covids)

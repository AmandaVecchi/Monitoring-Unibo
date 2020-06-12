#interpolation = to find measures not taken in the field
#download data (plot55) and put them n the LAB folder

setwd("C:/LAB/")
library(spatstat)

#import data (= table)
inp <- read.table("dati_plot55_LAST3.csv", sep=";", head=T)  #separate cloums by a semicolon + state first row is the header + name the table
head(inp)

#estimate of canopy cover where we do not have data
attack(inp) #state the data we are using
plot(X, Y)  #f we did not attach the set we should write: plot(inp$X, inp$Y)

#planar poit patter = expain what is the X and what is the Y, and the range of the two variables
summary(inp) #to see min and max value for the variables
inppp <- ppp(x=X, y=Y, c(716000,718000),c(4859000,4861000))

#marks = lable different point with the value of the canopy cover -> value we are using for our calculations. associate to each point the value of the canopy cover
# name(inp) to see names of columns
mark(inppp) <- Canopy.cov

#smooth = interpolation
canopy <- Smooth(inppp)

# plot the smooth and the points together
plot(canopy)
points(inpppp, col= "green")

marks(inppp) <- cop.lich.mean
lichs <- Smooth(inppp)
plot(lichs)
points(inppp)

#plot both plots together
par(mfrow=c(1,2))
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)

par(mfrow=c(1,3))
plot(canopy)
points(inppp)
plot(lichs)
points(inppp)
plot(Canopy.cov, cop.lich.mean, col="red", pch=19, cex=2)


####################
#working on another dataset

inp.psam <- read.table("dati_psammofile.csv", sep=";", head=T)
attach(inp.psam)
head(inp.psam)

dev.off() # close the par
plot(E, N) # plot the coordinates

summary(inp.psam)
inp.psam.ppp <- ppp(x=E,y=N,c(356450,372240),c(5059800,5064150))

names(inp.psam)
marks(inp.psam.ppp) <- C_org

C <- Smooth(inp.psam.ppp)
plot(C)
points(inp.psam.ppp)






























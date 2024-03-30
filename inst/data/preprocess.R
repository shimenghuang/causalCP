# Copyright (c) 2017
# All rights reserved.

library(zoo)

###
# Read and interpolate GDP data
###

# read in required dates from MPdata
MPdata <- read.csv("monetarypolicy.csv")
date <- as.yearmon(as.character(MPdata[,1]),format="%Y-%m")

# Swiss GDP
gdpch <- read.csv("GDPch_quarterly.csv")
gdpch[,1] <- as.yearmon(as.character(gdpch[,1]),format="%Y-%m-%d")

x <- as.numeric(gdpch[,1])
y <- gdpch[,2]
fit.gdpch <- splinefun(x,y,method="natural")
gdpch2 <- fit.gdpch(as.numeric(date))

plot(date,gdpch2,type="l")
points(x,gdpch[,2],col="red")

gdpch2 <- cbind(date,gdpch2)
write.csv(gdpch2,file="gdpch.csv")

# Euro GDP
gdpeu <- read.csv("GDPeu_quarterly.csv")
gdpeu[,1] <- as.yearmon(as.character(gdpeu[,1]),format="%Y-%m-%d")

x <- as.numeric(gdpeu[,1])
y <- gdpeu[,2]
fit.gdpeu <- splinefun(x,y,method="natural")
gdpeu2 <- fit.gdpeu(as.numeric(date))

plot(date,gdpeu2,type="l")
points(x,gdpeu[,2],col="red")

gdpeu2 <- cbind(date,gdpeu2)
write.csv(gdpeu2,file="gdpeu.csv")


###
# Read data (MPdata)
###

MPdata <- read.csv("monetarypolicy.csv")
MPdata[,1] <- as.yearmon(as.character(MPdata[,1]),format="%Y-%m")

save(MPdata,file="monetarypolicy.Rda")

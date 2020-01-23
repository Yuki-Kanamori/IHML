pred.y <- pred.a <- NULLgit      # true contribution
res.y <- res.a <- NULL       # contribution estimated by ML


# simulated data ------------------------------------------------
# SST
n.year <- 50     # number of year
year <- 1:n.year    # year index
a.sst <- 5    #   intercept of sst 
b.sst <- 0.4     #   trend of sst
sigma.sst <- 3     # variation in sst

sst <- rnorm(n.year, a.sst + b.sst*year, sigma.sst)     # sst realization
sst <- scale(sst)   # normalization


# SAL (Salinity)
n.area <- 50     # number of area
area <- 1:n.area    # area index
a.sal <- 10     #  intercept of sal
b.sal <- 0.2    #    trend of sal
sigma.sal <- 2    # variation in sal

sal <- rnorm(n.area, a.sal+b.sal*area, sigma.sal)       # sal realization
sal <- scale(sal)     # normalization



# Year from SST and SAL -----------------------------------------
intercept.year <- 2    # intercept of year effect
a.sst.to.year <- 0.2    # effect of sst for year
b.sst.to.year <- -0.3     # effect of sst^2 for year
a.sal.to.year <- 0.01      # effect of sal for year
b.sal.to.year <- -0.02      # effect of sal^2 for year
sigma.year <- 0.3     # additional variation in year effect

# mean year effect
mean.year.effect <- intercept.year + a.sst.to.year*sst + b.sst.to.year*sst^2 + a.sal.to.year*sal + b.sal.to.year*sal^2  

# mean year effect by only sst (marginalizing by SAL)
mean.year.effect.sst <- intercept.year+a.sst.to.year*sst+b.sst.to.year*sst^2+mean(a.sal.to.year*sal+b.sal.to.year*sal^2)

# mean year effect by only sal (marginalizing by sst)
mean.year.effect.sal <- intercept.year+mean(a.sst.to.year*sst+b.sst.to.year*sst^2)+a.sal.to.year*sal+b.sal.to.year*sal^2

# realized year effect = normalized mean year effect + variation
year.effect <- as.numeric(scale(mean.year.effect))+rnorm(n.year,0,sigma.year)    


# contribution of sst and sal for year
predict.year <- c(exp(-0.5*log(sum((year.effect-scale(mean.year.effect.sst))^2))),exp(-0.5*log(sum((year.effect-scale(mean.year.effect.sal))^2))))
# normalization
predict.year <- predict.year/sum(predict.year)



# Area from SST and SAL -----------------------------------------
intercept.area <- 2    # intercept of area effect
a.sst.to.area <- 0.02    # effect of sst for area
b.sst.to.area <- -0.03    # effect of sst^2 for area
a.sal.to.area <- 0.2    # effect of sal for area
b.sal.to.area <- 0.1      # effect of sal^2 for area
sigma.area <- 0.4     # additional variation in area effect

# mean area effect
mean.area.effect <- intercept.area+a.sst.to.area*sst+b.sst.to.area*sst^2+a.sal.to.area*sal+b.sal.to.area*sal^2

# mean area effect by only sst (marginalizing by sal)
mean.area.effect.sst <- intercept.area+a.sst.to.area*sst+b.sst.to.area*sst^2+mean(a.sal.to.area*sal+b.sal.to.area*sal^2)

# mean area effect by only sst (marginalizing by sal)
mean.area.effect.sal <- intercept.area+mean(a.sst.to.area*sst+b.sst.to.area*sst^2)+a.sal.to.area*sal+b.sal.to.area*sal^2

# realized area effect = normalized mean area effect + variation
area.effect <- as.numeric(scale(mean.area.effect))+rnorm(n.area,0,sigma.area) 


# contribution of sst and sal for area
predict.area <- c(exp(-0.5*log(sum((area.effect-scale(mean.area.effect.sst))^2))),exp(-0.5*log(sum((area.effect-scale(mean.area.effect.sal))^2))))
# normalization
predict.area <- predict.area/sum(predict.area)



# Year:Area -----------------------------------------------------
# variation in interaction
sigma.interact <- 0.25

# (normalized) year;area interaction with additional variation
interact <- as.numeric(scale(outer(year.effect, area.effect)))+matrix(rnorm(n.year*n.area,0,sigma.interact),nrow=n.year,ncol=n.area)    




# Density -------------------------------------------------------
# d(s,t)
dens <- matrix(year.effect,nrow=n.year,ncol=n.area)+matrix(area.effect,nrow=n.year,ncol=n.area,byrow=TRUE)+interact   # realized density

rownames(dens) <- 1:n.year
colnames(dens) <- 1:n.area


# d(*,t): year-specific density
dens.y <- rowMeans(dens)

# d(s,t)-d(*,t): area-specific density
dens.a <- sweep(dens,1,rowMeans(dens),FUN="-")



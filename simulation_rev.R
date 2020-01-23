## spatio-temporal model and machine learning
##

# Simulation

set.seed(1)

# Sim <- 100    # number of iterations
Sim <- 1

pred.y <- pred.a <- NULL     # true contribution
res.y <- res.a <- NULL       # contribution estimated by ML

for (i in 1:Sim){
  
  print(i)
  
  ## simulated data
  # SST  (Sea Surface Temperature)
  
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
  
  # Year from SST and SAL
  
  intercept.year <- 2    # intercept of year effect
  a.sst.to.year <- 0.2    # effect of sst for year
  b.sst.to.year <- -0.3     # effect of sst^2 for year
  a.sal.to.year <- 0.01      # effect of sal for year
  b.sal.to.year <- -0.02      # effect of sal^2 for year
  sigma.year <- 0.3     # additional variation in year effect
  
  mean.year.effect <- intercept.year+a.sst.to.year*sst+b.sst.to.year*sst^2+a.sal.to.year*sal+b.sal.to.year*sal^2    # mean year effect
  mean.year.effect.sst <- intercept.year+a.sst.to.year*sst+b.sst.to.year*sst^2+mean(a.sal.to.year*sal+b.sal.to.year*sal^2)   # mean year effect by only sst (marginalizing by sal)
  mean.year.effect.sal <- intercept.year+mean(a.sst.to.year*sst+b.sst.to.year*sst^2)+a.sal.to.year*sal+b.sal.to.year*sal^2    # mean year effect by only sal (marginalizing by sst)
  year.effect <- as.numeric(scale(mean.year.effect))+rnorm(n.year,0,sigma.year)    # realized year effect = normalized mean year effect + variation
  predict.year <- c(exp(-0.5*log(sum((year.effect-scale(mean.year.effect.sst))^2))),exp(-0.5*log(sum((year.effect-scale(mean.year.effect.sal))^2))))    # contribution of sst and sal for year
  predict.year <- predict.year/sum(predict.year)   # normalization
  
  # Area from SST and SAL
  
  intercept.area <- 2    # intercept of area effect
  a.sst.to.area <- 0.02    # effect of sst for area
  b.sst.to.area <- -0.03    # effect of sst^2 for area
  a.sal.to.area <- 0.2    # effect of sal for area
  b.sal.to.area <- 0.1      # effect of sal^2 for area
  sigma.area <- 0.4     # additional variation in area effect
  
  mean.area.effect <- intercept.area+a.sst.to.area*sst+b.sst.to.area*sst^2+a.sal.to.area*sal+b.sal.to.area*sal^2    # mean area effect
  mean.area.effect.sst <- intercept.area+a.sst.to.area*sst+b.sst.to.area*sst^2+mean(a.sal.to.area*sal+b.sal.to.area*sal^2)   # mean area effect by only sst (marginalizing by sal)
  mean.area.effect.sal <- intercept.area+mean(a.sst.to.area*sst+b.sst.to.area*sst^2)+a.sal.to.area*sal+b.sal.to.area*sal^2   # mean area effect by only sst (marginalizing by sal)
  area.effect <- as.numeric(scale(mean.area.effect))+rnorm(n.area,0,sigma.area)    # realized area effect = normalized mean area effect + variation
  predict.area <- c(exp(-0.5*log(sum((area.effect-scale(mean.area.effect.sst))^2))),exp(-0.5*log(sum((area.effect-scale(mean.area.effect.sal))^2))))    # contribution of sst and sal for area
  predict.area <- predict.area/sum(predict.area)   # normalization
  
  # Year:Area
  
  sigma.interact <- 0.25    # variation in interaction
  interact <- as.numeric(scale(outer(year.effect, area.effect)))+matrix(rnorm(n.year*n.area,0,sigma.interact),nrow=n.year,ncol=n.area)    # (normalized) year;area interaction with additional variation
  
  # density
  
  dens <- matrix(year.effect,nrow=n.year,ncol=n.area)+matrix(area.effect,nrow=n.year,ncol=n.area,byrow=TRUE)+interact   # realized density
  
  rownames(dens) <- 1:n.year
  colnames(dens) <- 1:n.area
  
  # year-specific and area-specific density
  
  dens.y <- rowMeans(dens)     # year-specific density 
  dens.a <- colMeans(sweep(dens,1,rowMeans(dens),FUN="-"))    # area-specific density
  
  # machine learning
  
  library(randomForest)    # activate random forest package
  
  mod.y1 <- randomForest(dens.y~sst)   # year-specific density predicted by sst
  mod.y2 <- randomForest(dens.y~sal)   # year-specific density predicted by sal
  mod.a1 <- randomForest(dens.a~sst)   # area-specific density predicted by sst
  mod.a2 <- randomForest(dens.a~sal)   # area-specific density predicted by sal
  
  contrib.y <- c(exp(-0.5*log(sum((dens.y-predict(mod.y1,type="response"))^2))),exp(-0.5*log(sum((dens.y-predict(mod.y2,type="response"))^2))))      # contribution of sst and sal for year estmated by random forest
  contrib.a <- c(exp(-0.5*log(sum((dens.a-predict(mod.a1,type="response"))^2))),exp(-0.5*log(sum((dens.a-predict(mod.a2,type="response"))^2))))      # contribution of sst and sal for area estmated by random forest
  names(contrib.y) <- names(contrib.a) <- c("SST","SAL")   # name of contribution
  
  contrib.y <- contrib.y/sum(contrib.y)    # normalization
  contrib.a <- contrib.a/sum(contrib.a)    # normalization
  
  # Outputs
  
  pred.y <- rbind(pred.y, predict.year)   # true contribution for year
  pred.a <- rbind(pred.a, predict.area)   # true contribution for area
  
  res.y <- rbind(res.y, contrib.y)    # estimated contribution for year
  res.a <- rbind(res.a, contrib.a)    # estimated contribution for area
}

colMeans(res.y)
colMeans(res.a)

year.bias <- (res.y-pred.y)/pred.y
area.bias <- (res.a-pred.a)/pred.a

res.bias <- data.frame(Year=year.bias[,1], Area=area.bias[,2])
boxplot(res.bias,col="brown")

dat1 <- data.frame(dens=scale(as.numeric(dens)),expand.grid(sst=sst,sal=sal))
test1 <- randomForest(dens~sst+sal,data=dat1)
test1
plot(dat1$dens,predict(test1,type="response"),xlab="Observed density",ylab="Predicted density",col="blue",cex=1.2)
abline(0,1,lwd=2,col="red",lty=2)
importance(test1)
varImpPlot(test1)
contrib.y
contrib.a
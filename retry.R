
# packages ------------------------------------------------------
library(CompRandFld)
library(RandomFields)
library(mapproj)
library(fields)
require(tidyverse)




pred.y <- pred.a <- NULL      # true contribution
res.y <- res.a <- NULL       # contribution estimated by ML


# simulated data ------------------------------------------------
# env(s,t) = env(t) + env(s) (+ env(s,t))
# 
# env(t) ~ N(mu, sigma)
# mu_sst = 5+0.4*Year / mu_sal = 10+0.2*Year
# 
# env(s) ~ MVN(0, Σ)
# ΣはMatérn関数(RFノイズはありで，localノイズは無し)

### SST(t)
n.year <- 50     # number of year
year <- 1:n.year    # year index
a.sst <- 5    #   intercept of sst 
b.sst <- 0.4     #   trend of sst
sigma.sst <- 3     # variation in sst

sst_t <- rnorm(n.year, a.sst + b.sst*year, sigma.sst)     # sst realization
sst_t = (sst_t - mean(sst_t))/sd(sst_t)

### SST(s)
# Define the spatial-coordinates of the points:
set.seed(261)
x = rep(1:15)
y = rep(1:15)

size = length(x)
size*size

# Simulation of a spatial Gaussian random field:
CorrelationParam("matern")
matern_sst <- RFsim(x, y, corrmodel = "matern", grid = TRUE, 
                    param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data

### heat map
mat_sst = matern_sst %>% data.frame() %>% gather(key = x, value = sst, 1:size) %>% 
  mutate(x = rep(1:size, each = size), y = rep(seq(1, size, 1), size))
summary(mat_sst)
mat_sst$sst = scale(mat_sst$sst)
colnames(mat_sst)

require(ggplot2)
g = ggplot(mat_sst, aes(x = x, y = y, fill = sst))
t = geom_tile()
c = scale_fill_gradientn(colours = c("black", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))
g+t+c+theme_bw()+labs(fill = "Scaled SST")

sst_s = mat_sst$sst

### SST(s,t)
# Define the spatial-coordinates of the points:
x = rep(1:15)
y = rep(1:15)

size = length(x)
size*size

# Simulation of a spatial Gaussian random field:
CorrelationParam("matern")
matern_sst <- RFsim(x, y, corrmodel = "matern", grid = TRUE, 
                    param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data

sst_st = matrix(NA, nrow = size*size, ncol = n.year)
rand = floor(runif(n.year, 0, 1000))
for(i in 1:n.year){
  set.seed(rand[i])
  matern = RFsim(x, y, corrmodel = "matern", grid = TRUE, 
              param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data
  eps = matern %>% data.frame() %>% gather(key = x, value = sst, 1:size) %>% 
    mutate(x = rep(1:size, each = size), y = rep(seq(1, size, 1), size))
  
  sst_st[, i] = eps$sst
}



### sal(t)
n.year <- 50     # number of year
year <- 1:n.year    # year index
a.sal <- 10     #  intercept of sal
b.sal <- 0.2    #    trend of sal
sigma.sal <- 2    # variation in sal

sal_t <- rnorm(n.year, a.sal+b.sal*year, sigma.sal)       # sal realization
sal_t = (sal_t - mean(sal_t))/sd(sal_t)   # normalization

# sal(s)
# Define the spatial-coordinates of the points:
set.seed(260)
x = rep(1:15)
y = rep(1:15)

size = length(x)
size*size

# Simulation of a spatial Gaussian random field:
CorrelationParam("matern")
matern_sal <- RFsim(x, y, corrmodel = "matern", grid = TRUE, 
                    param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data

### heat map
mat_sal = matern_sal %>% data.frame() %>% gather(key = x, value = sal, 1:size) %>% 
  mutate(x = rep(1:size, each = size), y = rep(seq(1, size, 1), size))
summary(mat_sal)
mat_sal$sal = scale(mat_sal$sal)
colnames(mat_sal)

require(ggplot2)
g = ggplot(mat_sal, aes(x = x, y = y, fill = sal))
t = geom_tile()
c = scale_fill_gradientn(colours = c("black", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))
g+t+c+theme_bw()+labs(fill = "Scaled sal")

sal_s = mat_sal$sal

### sal(s,t)
# Define the spatial-coordinates of the points:
x = rep(1:15)
y = rep(1:15)

size = length(x)
size*size

# Simulation of a spatial Gaussian random field:
CorrelationParam("matern")
matern_sst <- RFsim(x, y, corrmodel = "matern", grid = TRUE, 
                    param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data

sal_st = matrix(NA, nrow = size*size, ncol = n.year)
rand = floor(runif(n.year, 0, 1000))
for(i in 1:n.year){
  set.seed(rand[i])
  matern = RFsim(x, y, corrmodel = "matern", grid = TRUE, 
                 param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data
  eps = matern %>% data.frame() %>% gather(key = x, value = sal, 1:size) %>% 
    mutate(x = rep(1:size, each = size), y = rep(seq(1, size, 1), size))
  
  sal_st[, i] = eps$sal
}


# make SST and SAL ----------------------------------------------
# m_sst_t = matrix(sst_t, ncol = n.year, nrow = size*size, byrow = TRUE)
# m_sst_s = matrix(sst_s, ncol = n.year, nrow = size*size)
# summary(m_sst_s[,1] - sst_s)

A1 = 2
B1 = 0.5 #2だと0.643
C1 = 0.5
sst = A1*matrix(sst_t, ncol = n.year, nrow = size*size, byrow = TRUE) + B1*matrix(sst_s, ncol = n.year, nrow = size*size) + C1*sst_st
sst = (sst - mean(sst))/sd(sst)

A2 = 0.5
B2 = 0.5
C2 = 2
sal = A2*matrix(sal_t, ncol = n.year, nrow = size*size, byrow = TRUE) + B2*matrix(sal_s, ncol = n.year, nrow = size*size) + C2*sal_st
sal = (sal - mean(sal))/sd(sal)

# Year from SST and SAL -----------------------------------------
intercept.year <- 2    # intercept of year effect
a.sst.to.year <- 200    # effect of sst for year
b.sst.to.year <- -0.1     # effect of sst^2 for year
a.sal.to.year <- 0.001      # effect of sal for year
b.sal.to.year <- -0.002      # effect of sal^2 for year
sigma.year <- 0.3     # additional variation in year effect

# mean year effect
mean.year.effect <- intercept.year + a.sst.to.year*sst + b.sst.to.year*sst^2 + a.sal.to.year*sal + b.sal.to.year*sal^2  

# mean year effect by only sst (marginalizing by SAL)
mean.year.effect.sst <- intercept.year+a.sst.to.year*sst+b.sst.to.year*sst^2+mean(a.sal.to.year*sal+b.sal.to.year*sal^2)

# mean year effect by only sal (marginalizing by sst)
mean.year.effect.sal <- intercept.year+mean(a.sst.to.year*sst+b.sst.to.year*sst^2)+a.sal.to.year*sal+b.sal.to.year*sal^2

# realized year effect = normalized mean year effect + variation
year_noize = matrix(rnorm(size*size*n.year,0,sigma.year), ncol = n.year, nrow = size*size)
year.effect = (mean.year.effect - mean(mean.year.effect))/sd(mean.year.effect) + year_noize

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
area_noize = matrix(rnorm(size*size*n.year,0,sigma.year), ncol = n.year, nrow = size*size)
area.effect = (mean.area.effect - mean(mean.area.effect))/sd(mean.area.effect) + area_noize

# contribution of sst and sal for area (normalization)
predict.area = c(exp(-0.5*log(sum((area.effect-(mean.area.effect.sst - mean(mean.area.effect.sst))/sd(mean.area.effect.sst))^2))),exp(-0.5*log(sum((area.effect-(mean.area.effect.sal - mean(mean.area.effect.sal))/sd(mean.area.effect.sal))^2))))
predict.area <- predict.area/sum(predict.area)


# Density -------------------------------------------------------
# d(s,t): realized density
dens <- year.effect + area.effect

rownames(dens) <- 1:(size*size)
colnames(dens) <- 1:n.year


# d(*,t): year-specific density
dens.y <- colMeans(dens)

# d(s,t)-d(*,t): area-specific density
dens.a <- dens - dens.y



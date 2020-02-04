
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
#sst_t <- scale(sst_t)   # normalization

# SST(s)
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

summary(matern_sst)
sst_s = mat_sst %>% spread(key = x, value = sst) %>% select(-y) %>% as.matrix()

# sst(s,t)
sst_st = array(NA, dim = c(size,size,n.year))
for(i in 1:n.year){
  sst_st[,,i] = sst_t[i]*sst_s
}

# normalization
sst = (sst_st - mean(sst_st))/sd(sst_st)


# ####
# # ベクトル
# x <- c(1,2,3,4,5)
# zx <- (x - mean(x))/sd(x)
# zx
# scale(x)
# 
# # 行列
# m <- matrix(1:9,3,3)
# zm = (m-mean(m))/sd(m)
# zm
# 
# m2 = rep(1:9)
# zm2 = (m2-mean(m2))/sd(m2)
# zm2
# summary(zm2)
# 
# # 配列
# a = array(1:8, dim = c(2,2,2))
# za = (a-mean(a))/sd(a)
# za
# 
# a2 = rep(1:8)
# za2 = (a2-mean(a2))/sd(a2)
# za2
# summary(za2)
# ####



### sal(t)
n.year <- 50     # number of year
year <- 1:n.year    # year index
a.sal <- 10     #  intercept of sal
b.sal <- 0.2    #    trend of sal
sigma.sal <- 2    # variation in sal

sal_t <- rnorm(n.year, a.sal+b.sal*year, sigma.sal)       # sal realization
sal_t <- scale(sal_t)     # normalization

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

summary(matern_sal)
sal_s = mat_sal %>% spread(key = x, value = sal) %>% select(-y) %>% as.matrix()

# sal(s,t)
sal_st = array(NA, dim = c(size,size,n.year))
for(i in 1:n.year){
  sal_st[,,i] = sal_t[i]*sal_s
}

sal = (sal_st - mean(sal_st))/sd(sal_st)

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
# year.effect <- as.numeric(scale(mean.year.effect))+rnorm(n.year,0,sigma.year)    
# year_noize = rnorm(n.year,0,sigma.year) 
year_noize = array(rnorm(size*size*n.year,0,sigma.year), dim = c(size, size, n.year)) 
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
# area.effect <- as.numeric(scale(mean.area.effect))+rnorm(n.area,0,sigma.area) 
area_noize = array(rnorm(size*size*n.year, 0, sigma.area), dim = c(size, size, n.year))
area.effect = (mean.area.effect - mean(mean.area.effect))/sd(mean.area.effect) + area_noize

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
# d(s,t): realized density
dens <- matrix(year.effect,nrow=n.year,ncol=n.area)+matrix(area.effect,nrow=n.year,ncol=n.area,byrow=TRUE)+interact

rownames(dens) <- 1:n.year
colnames(dens) <- 1:n.area


# d(*,t): year-specific density
dens.y <- rowMeans(dens)

# d(s,t)-d(*,t): area-specific density
dens.a <- sweep(dens,1,rowMeans(dens),FUN="-")



# boosting ------------------------------------------------------
require(tidyverse)
require(xgboost)
require(Matrix)
require(caret)
#require(EIX)
#require(car)
require(mlr)
require(doParallel)
detectCores() # 4-cores
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)


# tune the hyper params -----------------------------------------
## not allow containing NA in data

# make data for tuning
d = data.frame(dens) %>% gather(key = year, value = dens, 1:50) %>% mutate(year = rep(1:50, each = 50), site = rep(seq(1,50,1), 50))

e1 = data.frame(sst) %>% mutate(year = rep(1:50))
e2 = data.frame(sal) %>% mutate(year = rep(1:50))

d = merge(d, e1, by = "year")
d = merge(d, e2, by = "year") 
d2 = d %>% select(-year, -site)

tuning = train(
  dens ~ .,
  data = d2,
  method = "xgbTree",
  preProcess = c("center", "scale"),
  trControl = trainControl(method = "cv"),
  tuneLength = 5
)

#results
best_tune = tuning$bestTune
plot(tuning)

params = list(
  booster           = 'gbtree',
  objective         = 'reg:linear',
  eval_metric       = 'mae',
  eta               = best_tune$eta,
  gamma             = best_tune$gamma,
  max_depth         = best_tune$max_depth,
  min_child_weight  = best_tune$min_child_weight,
  subsample         = best_tune$subsample,
  colsample_bytree  = best_tune$colsample_bytree
  )

summary(d2)
df = normalizeFeatures(d2, target = "dens")
summary(df)

full = xgboost(
  data = df %>% select(-dens) %>% as.matrix(),
  label = df$dens,
  nrounds = best_tune$nrounds,
  params = params)
pred_full = predict(model = full, data = df %>% select(-dens) %>% as.matrix())

  
model = get(paste0("model_sar_t", i))
data = df
data = data %>% filter(year == i) %>% select(-year,-month,-km,-knot_i)

pred = data.frame(predict(model, as.matrix(data[, -1])), data[,1])
colnames(pred) = c("pred_t", "obs")
pred = cbind(pred, filter(df, year == i) %>% select(year,month,km,knot_i)) %>% mutate(year = paste0(i))

pred_sar_t = rbind(pred_sar_t, pred)
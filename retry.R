
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

n.site = length(x)
n.site*n.site

# Simulation of a spatial Gaussian random field:
CorrelationParam("matern")
matern_sst <- RFsim(x, y, corrmodel = "matern", grid = TRUE, 
                    param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data

### heat map
mat_sst = matern_sst %>% data.frame() %>% gather(key = x, value = sst, 1:n.site) %>% 
  mutate(x = rep(1:n.site, each = n.site), y = rep(seq(1, n.site, 1), n.site))
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

n.site = length(x)
n.site*n.site

# Simulation of a spatial Gaussian random field:
CorrelationParam("matern")
matern_sst <- RFsim(x, y, corrmodel = "matern", grid = TRUE, 
                    param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data

sst_st = matrix(NA, nrow = n.site*n.site, ncol = n.year)
rand = floor(runif(n.year, 0, 1000))
for(i in 1:n.year){
  set.seed(rand[i])
  matern = RFsim(x, y, corrmodel = "matern", grid = TRUE, 
                 param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data
  eps = matern %>% data.frame() %>% gather(key = x, value = sst, 1:n.site) %>% 
    mutate(x = rep(1:n.site, each = n.site), y = rep(seq(1, n.site, 1), n.site))
  
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

n.site = length(x)
n.site*n.site

# Simulation of a spatial Gaussian random field:
CorrelationParam("matern")
matern_sal <- RFsim(x, y, corrmodel = "matern", grid = TRUE, 
                    param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data

### heat map
mat_sal = matern_sal %>% data.frame() %>% gather(key = x, value = sal, 1:n.site) %>% 
  mutate(x = rep(1:n.site, each = n.site), y = rep(seq(1, n.site, 1), n.site))
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

n.site = length(x)
n.site*n.site

# Simulation of a spatial Gaussian random field:
CorrelationParam("matern")
matern_sst <- RFsim(x, y, corrmodel = "matern", grid = TRUE, 
                    param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data

sal_st = matrix(NA, nrow = n.site*n.site, ncol = n.year)
rand = floor(runif(n.year, 0, 1000))
for(i in 1:n.year){
  set.seed(rand[i])
  matern = RFsim(x, y, corrmodel = "matern", grid = TRUE, 
                 param = list(smooth = 1, mean = 0, sill = 0.5, scale = 1, nugget = 0))$data
  eps = matern %>% data.frame() %>% gather(key = x, value = sal, 1:n.site) %>% 
    mutate(x = rep(1:n.site, each = n.site), y = rep(seq(1, n.site, 1), n.site))
  
  sal_st[, i] = eps$sal
}


# make SST and SAL ----------------------------------------------
# m_sst_t = matrix(sst_t, ncol = n.year, nrow = n.site*n.site, byrow = TRUE)
# m_sst_s = matrix(sst_s, ncol = n.year, nrow = n.site*n.site)
# summary(m_sst_s[,1] - sst_s)

A1 = 1
B1 = 1
C1 = 1
sst = A1*matrix(sst_t, ncol = n.year, nrow = n.site*n.site, byrow = TRUE) + B1*matrix(sst_s, ncol = n.year, nrow = n.site*n.site) + C1*sst_st
sst = (sst - mean(sst))/sd(sst)

A2 = 1
B2 = 1
C2 = 1
sal = A2*matrix(sal_t, ncol = n.year, nrow = n.site*n.site, byrow = TRUE) + B2*matrix(sal_s, ncol = n.year, nrow = n.site*n.site) + C2*sal_st
sal = (sal - mean(sal))/sd(sal)



# make Density -------------------------------------------------------
# d(s,t): realized density
gamma1 = 0.2
gamma2 = -0.3
gamma3 = 0.01
gamma4 = -0.02
sigma.year <- 0.3     # additional variation in year effect
sigma.area <- 0.4     # additional variation in area effect
  
sst_effect = gamma1*sst + gamma2*sst^2 + matrix(rnorm(n.site*n.site*n.year, 0, sigma.year), ncol = n.year, nrow = n.site*n.site)
sst_effect = (sst_effect - mean(sst_effect))/sd(sst_effect)
sal_effect = gamma3*sal + gamma4*sal^2 + matrix(rnorm(n.site*n.site*n.year, 0, sigma.area), ncol = n.year, nrow = n.site*n.site)
sal_effect = (sal_effect - mean(sal_effect))/sd(sal_effect)

b_sst_t = 0.3
b_sst_s = 0.01
b_sst_st = 0.01

b_sal_t = 0.02
b_sal_s = 0.5
b_sal_st = 0.5


# year effect terms
year_effect_sst = matrix(colMeans(sst_effect), ncol = n.year, nrow = n.site*n.site, byrow = TRUE)
year_effect_sal = matrix(colMeans(sal_effect), ncol = n.year, nrow = n.site*n.site, byrow = TRUE) 

# area effect terms
area_effect_sst = matrix(rowMeans(sst_effect), ncol = n.year, nrow = n.site*n.site)
area_effect_sal = matrix(rowMeans(sal_effect), ncol = n.year, nrow = n.site*n.site)

# interaction terms
int_effect_sst = (sst_effect - mean(sst_effect))
int_effect_sal = (sal_effect - mean(sal_effect))

# all
# dens = 
#   b_sst_t*(1/n.site*n.site)*(year_effect_sst - mean(year_effect_sst))/sd(year_effect_sst) + 
#   b_sst_s*(1/n.year)*(area_effect_sst - mean(area_effect_sst))/sd(area_effect_sst) + 
#   b_sst_st*(int_effect_sst - mean(int_effect_sst))/sd(int_effect_sst) +
#   b_sal_t*(1/n.site*n.site)*(year_effect_sal - mean(year_effect_sal))/sd(year_effect_sal) + 
#   b_sal_s*(1/n.year)*(area_effect_sal - mean(area_effect_sal))/sd(area_effect_sal) + 
#   b_sal_st*(int_effect_sal - mean(int_effect_sal))/sd(int_effect_sal)
dens = 
  b_sst_t*(year_effect_sst - mean(year_effect_sst))/sd(year_effect_sst) + 
  b_sst_s*(area_effect_sst - mean(area_effect_sst))/sd(area_effect_sst) + 
  b_sst_st*(int_effect_sst - mean(int_effect_sst))/sd(int_effect_sst) +
  b_sal_t*(year_effect_sal - mean(year_effect_sal))/sd(year_effect_sal) + 
  b_sal_s*(area_effect_sal - mean(area_effect_sal))/sd(area_effect_sal) + 
  b_sal_st*(int_effect_sal - mean(int_effect_sal))/sd(int_effect_sal)


# year effect -----------------------------------------
# mean year effect
mean.year.effect = dens #[225, 50]
# mean_effect = matrix(colMeans(sst_effect), ncol = n.year, nrow = n.site*n.site, byrow = TRUE)

# mean year effect by only sst (marginalizing by sal)
mean.year.effect.sst = 
  b_sst_t*(year_effect_sst - mean(year_effect_sst))/sd(year_effect_sst) + 
  mean(b_sst_s*(area_effect_sst - mean(area_effect_sst))/sd(area_effect_sst)) +
  mean(b_sst_st*(int_effect_sst - mean(int_effect_sst))/sd(int_effect_sst)) + 
  mean(b_sal_t*(year_effect_sal - mean(year_effect_sal))/sd(year_effect_sal)) + 
  mean(b_sal_s*(area_effect_sal - mean(area_effect_sal))/sd(area_effect_sal)) + 
  mean(b_sal_st*(int_effect_sal - mean(int_effect_sal))/sd(int_effect_sal))
s_y_sst = (mean.year.effect.sst - mean(mean.year.effect.sst))/sd(mean.year.effect.sst)

# mean year effect by only sal (marginalizing by sst)
mean.year.effect.sal = 
  mean(b_sst_t*(year_effect_sst - mean(year_effect_sst))/sd(year_effect_sst)) + 
  mean(b_sst_s*(area_effect_sst - mean(area_effect_sst))/sd(area_effect_sst)) + 
  mean(b_sst_st*(int_effect_sst - mean(int_effect_sst))/sd(int_effect_sst)) + 
  b_sal_t*(year_effect_sal - mean(year_effect_sal))/sd(year_effect_sal) + 
  mean(b_sal_s*(area_effect_sal - mean(area_effect_sal))/sd(area_effect_sal)) + 
  mean(b_sal_st*(int_effect_sal - mean(int_effect_sal))/sd(int_effect_sal))
s_y_sal = (mean.year.effect.sal - mean(mean.year.effect.sal))/sd(mean.year.effect.sal)

# year effect
year_noize = matrix(rnorm(n.site*n.site*n.year,0,sigma.year), ncol = n.year, nrow = n.site*n.site)
year.effect = (mean.year.effect - mean(mean.year.effect))/sd(mean.year.effect) + year_noize

# contribution of sst and sal for year
predict.year <- c(exp(-0.5*log(sum((year.effect-s_y_sst)^2))),exp(-0.5*log(sum((year.effect-s_y_sal)^2))))
# normalization
predict.year <- predict.year/sum(predict.year) #sst:sal = 0.55:0.45


# area effect -----------------------------------------
# mean area effect
mean.area.effect = dens #[225, 50]

# mean area effect by only sst (marginalizing by sal)
mean.area.effect.sst = 
  mean(b_sst_t*(year_effect_sst - mean(year_effect_sst))/sd(year_effect_sst)) + 
  b_sst_s*(area_effect_sst - mean(area_effect_sst))/sd(area_effect_sst) +
  b_sst_st*(int_effect_sst - mean(int_effect_sst))/sd(int_effect_sst) + 
  mean(b_sal_t*(year_effect_sal - mean(year_effect_sal))/sd(year_effect_sal)) + 
  mean(b_sal_s*(area_effect_sal - mean(area_effect_sal))/sd(area_effect_sal)) + 
  mean(b_sal_st*(int_effect_sal - mean(int_effect_sal))/sd(int_effect_sal))
s_a_sst = (mean.area.effect.sst - mean(mean.area.effect.sst))/sd(mean.area.effect.sst)

# mean year effect by only sal (marginalizing by sst)
mean.area.effect.sal = 
  mean(b_sst_t*(year_effect_sst - mean(year_effect_sst))/sd(year_effect_sst)) + 
  mean(b_sst_s*(area_effect_sst - mean(area_effect_sst))/sd(area_effect_sst)) + 
  mean(b_sst_st*(int_effect_sst - mean(int_effect_sst))/sd(int_effect_sst)) + 
  mean(b_sal_t*(year_effect_sal - mean(year_effect_sal))/sd(year_effect_sal)) + 
  b_sal_s*(area_effect_sal - mean(area_effect_sal))/sd(area_effect_sal) + 
  b_sal_st*(int_effect_sal - mean(int_effect_sal))/sd(int_effect_sal)
s_a_sal = (mean.area.effect.sal - mean(mean.area.effect.sal))/sd(mean.area.effect.sal)

# year effect
area_noize = matrix(rnorm(n.site*n.site*n.year,0,sigma.area), ncol = n.year, nrow = n.site*n.site)
area.effect = (mean.area.effect - mean(mean.area.effect))/sd(mean.area.effect) + area_noize

# contribution of sst and sal for year
predict.area <- c(exp(-0.5*log(sum((area.effect-s_a_sst)^2))),exp(-0.5*log(sum((area.effect-s_a_sal)^2))))
# normalization
predict.area <- predict.area/sum(predict.area) #sst:sal = 0.28:0.72



# Density -------------------------------------------------------
rownames(dens) = 1:(n.site*n.site)
colnames(dens) = 1:(n.year)
 
# year-specific and area-specific density
dens.y <- colMeans(dens)     #[1:50]
dens.a = matrix(NA, ncol = n.year, nrow = n.site*n.site) #[225,50]
for(i in 1:n.year){
  dens.a[, i] = dens[, i] - dens.y[i]
}
dens.y = matrix(dens.y, ncol = n.year, nrow = (n.site*n.site), byrow = TRUE) %>% data.frame()
dens.y = dens.y %>% gather(key = year, value = dens, 1:n.year) %>% mutate(year = rep(1:n.year, each = (n.site*n.site)), site = rep(seq(1, (n.site*n.site),1), n.year))
dens.a = dens.a %>% data.frame() %>% gather(key = year, value = dens, 1:n.year) %>% mutate(year = rep(1:n.year, each = (n.site*n.site)), site = rep(seq(1, (n.site*n.site),1), n.year))

# boosting ------------------------------------------------------
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
d = data.frame(dens) %>% gather(key = year, value = dens, 1:n.year) %>% mutate(year = rep(1:n.year, each = (n.site*n.site)), site = rep(seq(1, (n.site*n.site),1), n.year))

e1 = data.frame(sst) %>% gather(key = year, value = sst, 1:n.year) %>% mutate(year = rep(1:n.year, each = (n.site*n.site)), site = rep(seq(1, (n.site*n.site),1), n.year))
e2 = data.frame(sal) %>% gather(key = year, value = sal, 1:n.year) %>% mutate(year = rep(1:n.year, each = (n.site*n.site)), site = rep(seq(1, (n.site*n.site),1), n.year))

d = merge(d, e1, by = c("year", "site"))
d = merge(d, e2, by = c("year", "site")) 
d2 = d %>% select(-year, -site)
summary(d2)

set.seed(0)
tuning = caret::train(
  dens ~ ., 
  data = d2,
  method = "xgbTree", 
  preProcess = c('center', 'scale'),
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
sst = xgboost(
  data = df %>% select(-dens, -sal) %>% as.matrix(),
  label = df$dens,
  nrounds = best_tune$nrounds,
  params = params)
sal = xgboost(
  data = df %>% select(-dens, -sst) %>% as.matrix(),
  label = df$dens,
  nrounds = best_tune$nrounds,
  params = params)

# predict
model = full
data = df
pred_f = data.frame(predict(model, as.matrix(data[, -1])), data[,1])
colnames(pred_f) = c("pred", "obs")
plot(pred_f$obs, pred_f$pred)
cor(pred_f) #0.93

model = sst
data = df %>% select(-sal)
pred_sst = data.frame(predict(model, as.matrix(data[, -1])), data[,1])
colnames(pred_sst) = c("pred", "obs")
plot(pred_sst$obs, pred_sst$pred)
cor(pred_sst) #0.43

model = sal
data = df %>% select(-sst)
pred_sal = data.frame(predict(model, as.matrix(data[, -1])), data[,1])
colnames(pred_sal) = c("pred", "obs")
plot(pred_sal$obs, pred_sal$pred)
cor(pred_sal) #0.90



dens.y = merge(dens.y, e1, by = c("year", "site"))
dens.y = merge(dens.y, e2, by = c("year", "site")) 
dens.y = dens.y %>% select(-year, -site)

dens.a = merge(dens.a, e1, by = c("year", "site"))
dens.a = merge(dens.a, e2, by = c("year", "site")) 
dens.a = dens.a %>% select(-year, -site)

y_sst = xgboost(
  data = dens.y %>% select(-dens, -sal) %>% as.matrix(),
  label = dens.y$dens,
  nrounds = best_tune$nrounds,
  params = params)
y_sal = xgboost(
  data = dens.y %>% select(-dens, -sst) %>% as.matrix(),
  label = dens.y$dens,
  nrounds = best_tune$nrounds,
  params = params)

a_sst = xgboost(
  data = dens.a %>% select(-dens, -sal) %>% as.matrix(),
  label = dens.a$dens,
  nrounds = best_tune$nrounds,
  params = params)
a_sal = xgboost(
  data = dens.a %>% select(-dens, -sst) %>% as.matrix(),
  label = dens.a$dens,
  nrounds = best_tune$nrounds,
  params = params)

# predict
model = y_sst
data = dens.y %>% select(-sal)
pred_ysst = data.frame(predict(model, as.matrix(data[, -1])), data[,1])
colnames(pred_ysst) = c("pred", "obs")
plot(pred_ysst$obs, pred_ysst$pred)
cor(pred_ysst) #0.27

model = y_sal
data = dens.y %>% select(-sst)
pred_ysal = data.frame(predict(model, as.matrix(data[, -1])), data[,1])
colnames(pred_ysal) = c("pred", "obs")
plot(pred_ysal$obs, pred_ysal$pred)
cor(pred_ysal) #0.28

model = a_sst
data = dens.a %>% select(-sal)
pred_asst = data.frame(predict(model, as.matrix(data[, -1])), data[,1])
colnames(pred_asst) = c("pred", "obs")
plot(pred_asst$obs, pred_asst$pred)
cor(pred_asst) #0.30

model = a_sal
data = dens.a %>% select(-sst)
pred_asal = data.frame(predict(model, as.matrix(data[, -1])), data[,1])
colnames(pred_asal) = c("pred", "obs")
plot(pred_asal$obs, pred_asal$pred)
cor(pred_asal) #0.93


# contribution --------------------------------------------------
# contribution of sst and sal for full model
cont = c(exp(-0.5*log(sum((pred_f$obs - pred_sst$pred)^2))), exp(-0.5*log(sum((pred_f$obs - pred_sal$pred)^2))))
cont = cont/sum(cont) #sst:sal = 0.498:0.502

# contribution of sst and sal for year
cont_y = c(exp(-0.5*log(sum((pred_sst$obs - pred_ysst$pred)^2))), exp(-0.5*log(sum((pred_sal$obs - pred_ysal$pred)^2))))
cont_y = cont_y/sum(cont_y) #sst:sal = 0.497:0.503

# contribution of sst and sal for area
cont_a = c(exp(-0.5*log(sum((pred_sst$obs - pred_asst$pred)^2))), exp(-0.5*log(sum((pred_sal$obs - pred_asal$pred)^2))))
cont_a = cont_a/sum(cont_a) #sst:sal = 0.568:0.5432

names(cont_y) <- names(cont_a) <- c("SST","SAL")

# Outputs
pred.y <- rbind(pred.y, predict.year)   # true contribution for year
pred.a <- rbind(pred.a, predict.area)   # true contribution for area

res.y <- rbind(res.y, cont_y)    # estimated contribution for year
res.a <- rbind(res.a, cont_a)    # estimated contribution for area

colMeans(res.y)
colMeans(res.a)
predict.year
predict.area

year.bias <- (res.y-pred.y)/pred.y  
area.bias <- (res.a-pred.a)/pred.a 

res.bias <- data.frame(Year=year.bias[,1], Area=area.bias[,2])
boxplot(res.bias,col="brown")

# importance
Imp = xgb.importance(colnames(data), model = model)


# not revise
dat1 <- data.frame(dens=scale(as.numeric(dens)),expand.grid(sst=sst,sal=sal))
test1 <- randomForest(dens~sst+sal,data=dat1)
test1
plot(dat1$dens,predict(test1,type="response"),xlab="Observed density",ylab="Predicted density",col="blue",cex=1.2)
abline(0,1,lwd=2,col="red",lty=2)
importance(test1)
varImpPlot(test1)
contrib.y
contrib.a
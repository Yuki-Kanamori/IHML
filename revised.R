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
# => 空間パターンがモザイク状になることが気になる（生物の分布っぽくない）

# 変更
# env(s) = cos(s_1) - sin(s_2-2)

### SST(t)
set.seed(261)
n.year <- 50     # number of year
year <- 1:n.year    # year index
a.sst <- 5    #   intercept of sst 
b.sst <- 0.4     #   trend of sst
sigma.sst <- 3     # variation in sst

sst_t <- rnorm(n.year, a.sst + b.sst*year, sigma.sst)
sst_t = (sst_t - mean(sst_t))/sd(sst_t)

### SST(s)
# Define the spatial-coordinates of the points:
x = rep(1:15)
y = rep(1:15)
# n_pix = 300

n.site = length(x)
n.site*n.site


######################
setwd('/Users/Yuki/FRA/INLAren/spde-book-files')

## ----sett, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/prefsampl-'
)
library(splancs)

## ----mesh----------------------------------------------------------------
loc.d <- 3 * cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
# loc.d = cbind(rep(1:15, each = 15), rep(1:15))
plot(loc.d)
mesh <- inla.mesh.2d(loc.domain = loc.d, offset = c(0.3, 1), 
                     max.edge = c(0.3, 0.7), cutoff = 0.05)
plot(mesh)
mesh_loc = mesh$loc #492
nv <- mesh$n

## ------------------------------------------------------------------------
# sigma2x <- 0.2
# range <- 1.2
# nu <- 1
# library(spatstat)
# win = owin(c(0,3), c(0,3))
# 
# lg.s.c <- rLGCP('matern', im(1 * sst, xcol = x0,
#                              yrow = y0), var = sigma2x, scale = range / sqrt(8), 
#                 nu = nu)
######################

x0 = seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = n.site)
y0 = seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = n.site)
sst = outer(x0, y0, function(x, y) cos(1/2*x-2) - sin(1/2*y-1.2))
# lg.s.c <- rLGCP('matern', im(-0.5 * sst, xcol = x0,
#                              yrow = y0), var = sigma2x, scale = range / sqrt(8), nu = nu)
book.plot.field(list(x = x0, y = y0, z = sst))
summary(x0)
summary(y0)
# book.plot.field(list(x = x0, y = y0, z = sst), xlim = c(0, 3), ylim = c(min(y0), 4))

sst = sst %>% data.frame() %>% mutate(x0 = x0)　%>% gather(key = x, value = sst, 1:n.site)
tag = data.frame(x = paste0("X", rep(1:n.site)), y0 = y0)
sst = left_join(sst, tag, by = "x")
# sst$sst = sst$sst + rnorm(n.site*n.site, 0, sigma.sst)
sst$sst = scale(sst$sst)
colnames(sst)

require(ggplot2)
g = ggplot(sst, aes(x = x0, y = y0, fill = sst))
t = geom_tile()
c = scale_fill_gradientn(colours = c("black", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))
g+t+c+theme_bw()+labs(fill = "SST")

sst_s = sst$sst




### sal(t)
set.seed(261)
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

loc.d <- 3 * cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
# loc.d = cbind(rep(1:15, each = 15), rep(1:15))
plot(loc.d)
mesh <- inla.mesh.2d(loc.domain = loc.d, offset = c(0.3, 1), 
                     max.edge = c(0.3, 0.7), cutoff = 0.05)
plot(mesh)
mesh_loc = mesh$loc #492
nv <- mesh$n


x0 = seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = n.site)
y0 = seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = n.site)
sal = outer(x0, y0, function(x, y) cos(1/2*x) - sin(1/2*y-1.2))
# lg.s.c <- rLGCP('matern', im(-0.5 * sst, xcol = x0,
#                              yrow = y0), var = sigma2x, scale = range / sqrt(8), nu = nu)
book.plot.field(list(x = x0, y = y0, z = sal))
summary(x0)
summary(y0)
# book.plot.field(list(x = x0, y = y0, z = sst), xlim = c(0, 3), ylim = c(min(y0), 4))

sal = sal %>% data.frame() %>% mutate(x0 = x0)　%>% gather(key = x, value = sal, 1:n.site)
tag = data.frame(x = paste0("X", rep(1:n.site)), y0 = y0)
sal = left_join(sal, tag, by = "x")
# sst$sst = sst$sst + rnorm(n.site*n.site, 0, sigma.sst)
sal$sal = scale(sal$sal)
colnames(sal)

require(ggplot2)
g = ggplot(sal, aes(x = x0, y = y0, fill = sal))
t = geom_tile()
c = scale_fill_gradientn(colours = c("black", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))
g+t+c+theme_bw()+labs(fill = "SST")

sal_s = sal$sal


# make SST and SAL ----------------------------------------------
# m_sst_t = matrix(sst_t, ncol = n.year, nrow = n.site*n.site, byrow = TRUE)
# m_sst_s = matrix(sst_s, ncol = n.year, nrow = n.site*n.site)
# summary(m_sst_s[,1] - sst_s)

A1 = 1
B1 = 1
C1 = 1
sst = A1*matrix(sst_t, ncol = n.year, nrow = n.site*n.site, byrow = TRUE) + B1*matrix(sst_s, ncol = n.year, nrow = n.site*n.site) 
# + C1*sst_st
sst = (sst - mean(sst))/sd(sst)

A2 = 1
B2 = 1
C2 = 1
sal = A2*matrix(sal_t, ncol = n.year, nrow = n.site*n.site, byrow = TRUE) + B2*matrix(sal_s, ncol = n.year, nrow = n.site*n.site)
# + C2*sal_st
sal = (sal - mean(sal))/sd(sal)


# make Density -------------------------------------------------------
# d(s,t): realized density
gamma1 = 0.2
gamma2 = -0.3
gamma3 = 0.01
gamma4 = -0.02

sst_effect = gamma1*sst + gamma2*sst^2
sst_effect = (sst_effect - mean(sst_effect))/sd(sst_effect)
sal_effect = gamma3*sal + gamma4*sal^2
sal_effect = (sal_effect - mean(sal_effect))/sd(sal_effect)

b_sst_t = 0.3
b_sst_s = 0.01
b_sst_st = 0.01

b_sal_t = 0.02
b_sal_s = 0.5
b_sal_st = 0.5

sigma.year <- 0.3     # additional variation in year effect
sigma.area <- 0.4     # additional variation in area effect

# year effect terms
year_effect_sst = matrix(colMeans(sst_effect + matrix(rnorm(n.site*n.site*n.year, 0, sigma.year), ncol = n.year, nrow = n.site*n.site)), ncol = n.year, nrow = n.site*n.site, byrow = TRUE) 
year_effect_sal = matrix(colMeans(sal_effect + matrix(rnorm(n.site*n.site*n.year, 0, sigma.year), ncol = n.year, nrow = n.site*n.site)), ncol = n.year, nrow = n.site*n.site, byrow = TRUE) 

# area effect terms
area_effect_sst = matrix(b_sst_s*(1/n.year)*rowMeans(sst_effect + matrix(rnorm(n.site*n.site*n.year, 0, sigma.area), ncol = n.year, nrow = n.site*n.site)), ncol = n.year, nrow = n.site*n.site)
area_effect_sal = matrix(b_sal_s*(1/n.year)*rowMeans(sal_effect + matrix(rnorm(n.site*n.site*n.year, 0, sigma.area), ncol = n.year, nrow = n.site*n.site)), ncol = n.year, nrow = n.site*n.site)

# interaction terms
int_effect_sst = (sst_effect - mean(sst_effect))
int_effect_sal = (sal_effect - mean(sal_effect))

# all
dens = 
  b_sst_t*(1/n.site*n.site)*(year_effect_sst - mean(year_effect_sst))/sd(year_effect_sst) + 
  b_sst_s*(1/n.year)*(area_effect_sst - mean(area_effect_sst))/sd(area_effect_sst) + 
  b_sst_st*(int_effect_sst - mean(int_effect_sst))/sd(int_effect_sst) +
  b_sal_t*(1/n.site*n.site)*(year_effect_sal - mean(year_effect_sal))/sd(year_effect_sal) + 
  b_sal_s*(1/n.year)*(area_effect_sal - mean(area_effect_sal))/sd(area_effect_sal) + 
  b_sal_st*(int_effect_sal - mean(int_effect_sal))/sd(int_effect_sal)


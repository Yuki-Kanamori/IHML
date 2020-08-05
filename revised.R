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
n.year <- 50     # number of year
year <- 1:n.year    # year index
a.sst <- 5    #   intercept of sst 
b.sst <- 0.4     #   trend of sst
sigma.sst <- 3     # variation in sst

sst_t <- rnorm(n.year, a.sst + b.sst*year, sigma.sst)
sst_t = (sst_t - mean(sst_t))/sd(sst_t)

### SST(s)
# Define the spatial-coordinates of the points:
set.seed(1)
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

sst = sst %>% data.frame() %>% mutate(x0 = x0) %>% gather(key = x, value = sst, 1:n.site-1)
tag = data.frame(x = paste0("X", rep(1:n.site)), y0 = y0)
sst = left_join(sst, tag, by = "x")
sst$sst = scale(sst$sst)
colnames(sst)

require(ggplot2)
g = ggplot(sst, aes(x = x0, y = y0, fill = sst))
t = geom_tile()
c = scale_fill_gradientn(colours = c("black", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))
g+t+c+theme_bw()+labs(fill = "SST")

sst_s = sst$sst























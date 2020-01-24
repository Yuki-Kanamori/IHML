# NOT RUN {
library(CompRandFld)
library(RandomFields)
library(mapproj)
library(fields)
require(tidyverse)

################################################################
###
### Example 1. Simulation of a Gaussian random field.
### Gaussian random fields with Whittle-Matern correlation.
### One spatial replication.
###
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(500, 0, 2)
y <- runif(500, 0, 2)

set.seed(261)
# Simulation of a spatial Gaussian random field:
data <- RFsim(x, y, corrmodel="matern", param=list(smooth=0.5,
                                                   mean=0,sill=1,scale=0.2,nugget=0))$data



################################################################
###
### Example 4. Simulation of a Gaussian random field.
### with double exponential correlation.
### One spatio-temporal replication.
###
###
###############################################################

# Define the spatial-coordinates of the points:
x <- seq(0, 3, 0.1)
y <- seq(0, 3, 0.1)
# Define the temporal-coordinates:
times <- seq(1, 10, 1)
#
# Simulation of a spatial Gaussian random field:
sim <- RFsim(x, y, times, corrmodel = "exp_exp", grid = TRUE,
             param=list(nugget = 0, mean = 0, scale_s = 0.3,
                        scale_t = 0.5, sill = 1))$data
# Spatial simulated data at first temporal instant
sim[,,1]
dim(sim)[3]

df = NULL
for(i in 1:dim(sim)[3]){
  data = sim[,,i] %>% data.frame() %>% gather(key = x, value = logdens, 1:31) %>% 
    mutate(x = rep(1:31, each = 31), y = rep(seq(1, 31, 1), 31), year = paste0(i));
  
  df = rbind(df, data)
}

mode(df$year)
df$year = as.numeric(df$year)
#df$year = factor(df$year, levels = c("1","2","3","4","5","6","7","8","9","10"))

require(ggplot2)
g = ggplot(df, aes(x = x, y = y, fill = logdens))
t = geom_tile()
f = facet_wrap(~ year, ncol = 5)
c = scale_fill_gradientn(colours = c("black", "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))
g+t+c+f+theme_bw()

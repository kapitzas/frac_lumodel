rm(list = ls())

require(raster)
require(extraDistr)
source("./R/functions.R")

#---------------------------#
#### 1. DATA PREPERATION ####
#---------------------------#

# 1.a Load data

data_path <- file.path(getwd(), "data", "data_ama")
dat <- readRDS(file.path(data_path, "cov.rds")) #dynamic bioclimatic variables
lu_all <- readRDS(file.path(data_path, "lu.rds"))
mask <- readRDS(file.path(data_path, "mask_ama.rds")) #country mask

ts <- 1991 + c(1, 5, 10, 15, 20, 25, 27)
K <- 13
n <- nrow(lu_all)
lu <- lu_all[,1:K]

# 2.b Neighbourhood
weights <- list(matrix(1/9, 3, 3, byrow= TRUE)) #size of window
weights <- rep(weights, length.out = K)
ln <- neighbourhood(lu, c(1:K), weights, mask, enr = TRUE)

ln <- scale(ln)

cent <- attr(ln,c("scaled:center"))
scal <- attr(ln,c("scaled:scale"))

# 2.c Correlation analysis and susbetting
dat <- scale(dat)
data <- cbind(dat, ln)

preds <- colnames(correlations(data, sub = nrow(data)))
data <- data[,colnames(data)%in%preds]

#----------------------------#
#### 2. SUITABILITY MODEL ####
#----------------------------#

form <- paste(preds, collapse = "+")
subs_mod <- sample(1:nrow(lu), 20000)
suitmod <- suitmodel(form = form, lu = lu[subs_mod,], data = data[subs_mod,], resolution = 10000, model = FALSE, maxit = 10000)
sm <- predict(suitmod, newdata = data, type = "probs")

#-----------------#
#### 3. DEMAND ####
#-----------------#

demands <- demand(landuse = lu_all, ts = ts, inds = NULL, k = K, type = "mean")[,1:(K+1)]

#---------------------#
#### 4. SIMULATIONS####
#---------------------#

# Simulation parameters
params <- list(
  stepsi = 1,
  max_dev = 1,
  resolution = 1000000,
  max_it = 20000,
  ch_thresh = c(0.5, 0.5, #crop, crop_natveg_mosaic, 
                0.5, 0.5, 0.5, #tree_co, tree_mixed, treeshrub_herb_mosaic, 
                0.5, 0.5,  0.5, #shrub, grass, spares, 
                0.5, 0.5, #tree_water, shrub_water, 
                0.5, #urban, 
                0.5, 0.5),#bare, snow_ice
  no_change = 13
)

# 4.a Fully Naive model simulation

#Demand
dmd <- demands
dmd <- dmd[which(dmd[,1]%in%ts),-1]

#Prepare variables
lu_out <- lu
lu_ts <- list()
lu_pred <- matrix(NA, nrow(lu), ncol(lu))
colnames(lu_pred) <- colnames(lu)

#Simulate time series
for(i in 1:(length(ts)-1)){
  cat(paste0('\n', "Predicting suitability for time step ", i))
  
  #Random suitability
  sm <- matrix(runif(K*nrow(lu)), ncol = K, nrow = nrow(lu))
  sm <- sm/rowSums(sm)
  
  #Time step demand
  dmd_t0 <- dmd[i,]
  if(i > 1){
    dmd_t0 <- colMeans(lu_ts[[i-1]])
  }
  dmd_ts <- rbind(dmd_t0, dmd[i+1,])
  
  #Allocation
  lu_pred <- allocation(lu = lu_out,
                        ln = ln,
                        sm = sm, 
                        params = params, 
                        dmd = dmd_ts,
                        constraint = NULL)
  
  cat('\n')
  lu_ts[[i]] <- lu_out <- lu_pred
}

saveRDS(lu_ts, file = file.path("outputs", "preds_naive.rds"))

# 4.b semi-naive model

#Demand
dmd <- demands
dmd <- dmd[which(dmd[,1]%in%ts),-1]

#Prepare variables
lu_out <- lu
lu_ts <- list()
lu_pred <- matrix(NA, nrow(lu), ncol(lu))
colnames(lu_pred) <- colnames(lu)

#Simulate time series
for(i in 1:(length(ts)-1)){
  cat(paste0('\n', "Predicting suitability for time step ", i))
  
  #Predict suitability model
  ln <- neighbourhood(lu_out, 1:K, weights, mask, enr = TRUE)
  sm <- predict(suitmod, newdata = cbind(dat, scale(ln, center = cent, scale = scal)), type = "probs")
  
  #Time step demand
  dmd_t0 <- dmd[i,]
  if(i > 1){
    dmd_t0 <- colMeans(lu_ts[[i-1]])
  }
  dmd_ts <- rbind(dmd_t0, dmd[i+1,])
  
  #Allocation
  lu_pred <- allocation(lu = lu_out,
                        ln = ln,
                        sm = sm, 
                        params = params, 
                        dmd = dmd_ts,
                        constraint = NULL)
  
  cat('\n')
  lu_ts[[i]] <- lu_out <- lu_pred
}

saveRDS(lu_ts, file = file.path("outputs", "preds_semi.rds"))

#4. c Full model

#Demand
dmd <- demands
dmd <- dmd[which(dmd[,1]%in%ts),-1]

#Prepare variables
lu_out <- lu
lu_ts <- list()
lu_pred <- matrix(NA, nrow(lu), ncol(lu))
colnames(lu_pred) <- colnames(lu)

#Simulate time series
for(i in 1:(length(ts)-1)){
  cat(paste0('\n', "Predicting suitability for time step ", i))
  
  #Predict suitability model
  ln <- neighbourhood(lu_out, 1:K, weights, mask, enr = TRUE)
  #sm <- predict(suitmod, newdata = cbind(dat, scale(ln, center = cent, scale = scal)), type = "probs")
  sm <- matrix(runif(K*nrow(lu)), ncol = K, nrow = nrow(lu))
  sm <- sm/rowSums(sm)
  
  #Time step demand
  dmd_t0 <- dmd[i,]
  if(i > 1){
    dmd_t0 <- colMeans(lu_ts[[i-1]])
  }
  dmd_ts <- rbind(dmd_t0, dmd[i+1,])
  
  
  #Allocation
  lu_pred <- allocation(lu = lu_out,
                        ln = ln,
                        sm = sm, 
                        params = params, 
                        dmd = dmd_ts,
                        constraint = "all")
  
  cat('\n')
  lu_ts[[i]] <- lu_out <- lu_pred
}

saveRDS(lu_ts, file = file.path("outputs", "preds_full.rds"))





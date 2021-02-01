rm(list = ls())
require(raster)
require(extraDistr)
require(wrswoR)
library(flutes)

#---------------------------#
#### 1. DATA PREPERATION ####
#---------------------------#

# 1.a Load data (outputs of preprocessing script)

data_path <- file.path(getwd(), "data", "data_ama")
dat <- readRDS(file.path(data_path, "cov.rds")) #dynamic bioclimatic variables head(dat)
lu_all <- readRDS(file.path(data_path, "lu.rds"))
mask <- readRDS(file.path(data_path, "mask_ama.rds")) #country mask
inds <- which(!is.na(mask[]))
colnames(dat)

# 1.b define variables to prepare data for validaiton

# time steps
ts <- 1991 + c(1, 5, 10, 15, 20, 25, 27)
ts_inds <- 1:length(ts)

# number of lu classes
K <- ncol(lu_all)/length(ts)

# number of cells
n <- nrow(lu_all)

# reorganise observed land use data into a list of data frames (each list element contains one validation time step, columns are land uses, rows are cells)
lu_obs <- list()
start <- seq(1,ncol(lu_all), by = K)
for(i in 1:length(start)){
  lu_obs[[i]] <- as.data.frame(lu_all[,start[i]:(start[i]+(K-1))])
  names(lu_obs[[i]]) <- paste0("lu", 1:K, ".obs")
  print(i)
}

# first observed land use time step to parametrise the model
lu <- lu_all[,1:K]

# 1.c Calculate neighbourhood covariate (as "enrichment factors").
weights <- list(matrix(1/9, 3, 3, byrow= TRUE))
weights <- rep(weights, length.out = K)

# function flutes::neighbourhood described in flutes GitHub package
ln <- neighbourhood(lu, c(1:K), weights, mask, enr = TRUE)

# scale
ln <- scale(ln)

# store scaling parameters for scaling neighbourhood covariate in future time steps

cent <- attr(ln,c("scaled:center"))
scal <- attr(ln,c("scaled:scale"))

# 1.d Correlation analysis and susbetting

# exclude protected areas from correlation analysis
pa <- dat[,which(grepl("PA", colnames(dat)))]
pa <- pa[,1]
dat <- dat[,-which(grepl("PA", colnames(dat)))]

# scale covariates
dat <- scale(dat)
data <- cbind(dat, ln)

# estimate correlation coefficients and reduce predictor set  (function "flutes::correlations" part of flutes GitHub package)
preds <- colnames(correlations(data, sub = nrow(data)))
data <- data[,colnames(data)%in%preds]

# extract remaining predictor names for the model forumla below
preds <- colnames(data)


# 1. e extract how many cells changed from 0 to containing a land use (this is the data in Table 4)
ch_ma <- matrix(NA, nrow = 6, ncol = K)
for(i in 2:7){
  for(j in 1:K){
    ch_ma[i-1, j] <- length(lu_obs[[i-1]][which(lu_obs[[i-1]][,j] == 0 & lu_obs[[i]][,j] != 0),i])/n * 100
  }
}

lu_classes <-  c("Cro", "CrM", "For", "Gra", "Shr", "Wet", "Urb", "Oth", "Wat")
colnames(ch_ma) <- lu_classes
rownames(ch_ma) <- ts[-1]
ch_ma <- t(ch_ma)
ch_ma <- cbind(ch_ma, "mean" = rowMeans(ch_ma))
saveRDS(ch_ma, file = file.path("outputs", "lu_newestablishment.rds"))

#1. f demand changes
demands <- demand(landuse = lu_all, ts = ts, k = K, type = "mean")[,1:(K+1)]
demands[,1+K] <- rep(demands[1,1+K], nrow(demands))
demands[,-1] <- demands[,-1]/rowSums(demands[,-1])
saveRDS(demands, file = file.path(data_path, "demands.rds"))

#----------------------------#
#### 2. SUITABILITY MODEL ####
#----------------------------#

# 2. a build model formula
form <- paste(preds, collapse = "+")

# 2. b sample random subset of cells to build model on
subs_mod <- sample(1:nrow(lu), 50000)

# 2. c build model (function flutes::suitmodel in GitHub package)
suitmod <- suitmodel(form = form, lu = lu[subs_mod,], data = data[subs_mod,], resolution = 10000, model = FALSE, maxit = 10000, decay = 0.01)

#2. d predict suitability surfaces for current time step
sm <- predict(suitmod, newdata = data, type = "probs")

#---------------------#
#### 3. SIMULATIONS####
#---------------------#
# We simulate under three parametrization levels: naive, semi-naive, full to assess the respective contributions of each parametrization stage to predictive performance

# 3. a Define Simulation parameters
# max dev is the % deviation we allow between allocated supply and prescribed demand, resolution is the integer count to which fractions are converted for multinomial model. i.e. a fraction of 0.2 would become round(0.2 * resolution), growth is the % of the landscape (in terms of cells) where land use was 0 in a land use and contained a fraction of that land-use in the next time step, averaged thgroygh the observed time series. no_change is a vector of land use types that aren't allowed to change at all. Here we only do that for one land use "water", assuming that fresh water coverage in the amazon stays approximately the same.

params <- list(
  max_dev = 1,
  resolution = 1000000,
  growth = ch_ma[,7],
  no_change = 9
)

# 3.b Fully Naive model simulation

# Demand
dmd <- demands
dmd <- dmd[which(dmd[,1]%in%ts),-1]

# Prepare variables
lu_out <- lu
lu_ts <- list()
lu_suit <- list()
lu_pred <- matrix(NA, nrow(lu), ncol(lu))
colnames(lu_pred) <- colnames(lu)

# Simulate time series
for(i in 1:(length(ts)-1)){
  cat(paste0('\n', "Predicting suitability for time step ", i))
  
  # Random suitability ("naive information")
  sm <- matrix(runif(K*nrow(lu)), ncol = K, nrow = nrow(lu))
  sm <- sm/rowSums(sm)
  
  # Time step demand
  dmd_t0 <- dmd[i,]
  if(i > 1){
    dmd_t0 <- colMeans(lu_ts[[i-1]])
  }
  dmd_ts <- rbind(dmd_t0, dmd[i+1,])
  
  # Allocation (naive, without suitability maps or constraints)
  lu_pred <- allocation(lu = lu_out,
                        ln = ln,
                        sm = sm, 
                        params = params, 
                        dmd = dmd_ts,
                        constraint = FALSE,
                        pa = NULL)
  
  cat('\n')
  lu_ts[[i]] <- lu_out <- lu_pred
  lu_suit[[i]] <- sm
}

# check that all predictions are within defined precision limits
all(abs((do.call("rbind", lapply(lu_ts, FUN = colMeans)) - dmd[-1,])/dmd[-1,] * 100) < 1) #Must be TRUE
inds_pa <- which(pa == 0)
all(unlist(lapply(lu_ts, FUN = function(x) {all(x[inds_pa,]== lu[inds_pa,])}))) #must be FALSE if pa = NULL

saveRDS(lu_ts, file = file.path("outputs", "preds_naive.rds"))
saveRDS(lu_suit, file = file.path("outputs", "suit_naive.rds"))

# 3.c semi-naive model

# Demand
dmd <- demands
dmd <- dmd[which(dmd[,1]%in%ts),-1]

# Prepare variables
lu_out <- lu
lu_ts <- list()
lu_suit <- list()
lu_pred <- matrix(NA, nrow(lu), ncol(lu))
colnames(lu_pred) <- colnames(lu)

# Simulate time series
for(i in 1:(length(ts)-1)){
  cat(paste0('\n', "Predicting suitability for time step ", i))
  
  # Predict suitability model
  ln <- neighbourhood(lu_out, 1:K, weights, mask, enr = TRUE)
  sm <- predict(suitmod, newdata = cbind(dat, scale(ln, center = cent, scale = scal)), type = "probs")
  
  # Time step demand
  dmd_t0 <- dmd[i,]
  if(i > 1){
    dmd_t0 <- colMeans(lu_ts[[i-1]])
  }
  dmd_ts <- rbind(dmd_t0, dmd[i+1,])
  
  # Allocation (semi-naive, with proper suitability maps, but without any other constraints)
  lu_pred <- allocation(lu = lu_out,
                        ln = ln,
                        sm = sm, 
                        params = params, 
                        dmd = dmd_ts,
                        constraint = FALSE,
                        pa = NULL)
  
  cat('\n')
  lu_ts[[i]] <- lu_out <- lu_pred
  lu_suit[[i]] <- sm
}

# check that all predictions are within defined precision limits
all(abs((do.call("rbind", lapply(lu_ts, FUN = colMeans)) - dmd[-1,])/dmd[-1,] * 100) < 1) #Must be TRUE
inds_pa <- which(pa == 0)
all(unlist(lapply(lu_ts, FUN = function(x) {all(x[inds_pa,]== lu[inds_pa,])}))) #must be FALSE if pa = NULL

saveRDS(lu_ts, file = file.path("outputs", "preds_semi.rds"))
saveRDS(lu_suit, file = file.path("outputs", "suit_semi.rds"))

#3. d Full model

# Demand
dmd <- demands
dmd <- dmd[which(dmd[,1]%in%ts),-1]

# Prepare variables
lu_out <- lu
lu_ts <- list()
lu_suit <- list()
lu_pred <- matrix(NA, nrow(lu), ncol(lu))
colnames(lu_pred) <- colnames(lu)

# Simulate time series
for(i in 1:(length(ts)-1)){
  cat(paste0('\n', "Predicting suitability for time step ", i))
  
  # Predict suitability model
  ln <- neighbourhood(lu_out, 1:K, weights, mask, enr = TRUE)
  sm <- predict(suitmod, newdata = cbind(dat, scale(ln, center = cent, scale = scal)), type = "probs")
  
  # Time step demand
  dmd_t0 <- dmd[i,]
  if(i > 1){
    dmd_t0 <- colMeans(lu_ts[[i-1]])
  }
  dmd_ts <- rbind(dmd_t0, dmd[i+1,])
  
  
  # Allocation (full, with suitability maps and constraints)
  lu_pred <- allocation(lu = lu_out,
                        ln = ln,
                        sm = sm, 
                        params = params, 
                        dmd = dmd_ts,
                        constraint = TRUE,
                        pa = pa)
  
  cat('\n')
  lu_ts[[i]] <- lu_out <- lu_pred
  lu_suit[[i]] <- sm
}

# check that all predictions are within defined precision limits
all(abs((do.call("rbind", lapply(lu_ts, FUN = colMeans)) - dmd[-1,])/dmd[-1,] * 100) < 1) #Must be TRUE
inds_pa <- which(pa == 0)
all(unlist(lapply(lu_ts, FUN = function(x) {all(x[inds_pa,]== lu[inds_pa,])}))) #must be TRUE if pa = NULL

saveRDS(lu_ts, file = file.path("outputs", "preds_full.rds"))
saveRDS(lu_suit, file = file.path("outputs", "suit_full.rds"))
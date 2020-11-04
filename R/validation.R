rm(list = ls())

require(raster)
require(extraDistr)
require(wrswoR)
devtools::install("/Users/simon/OneDrive - The University of Melbourne/PhD/packages/flutes")
library(flutes)
suitmodel
#---------------------------#
#### 1. DATA PREPERATION ####
#---------------------------#

# 1.a Load data

data_path <- file.path(getwd(), "data", "data_ama")
dat <- readRDS(file.path(data_path, "cov.rds")) #dynamic bioclimatic variables head(dat)
lu_all <- readRDS(file.path(data_path, "lu.rds"))
mask <- readRDS(file.path(data_path, "mask_ama.rds")) #country mask
inds <- which(!is.na(mask[]))
colnames(dat)

ts <- 1991 + c(1, 5, 10, 15, 20, 25, 27)
K <- ncol(lu_all)/length(ts)
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
pa <- dat[,which(grepl("PA", colnames(dat)))]
pa <- pa[,1]
dat <- dat[,-which(grepl("PA", colnames(dat)))]
dat <- scale(dat)
data <- cbind(dat, ln)

preds <- colnames(correlations(data, sub = nrow(data)))
data <- data[,colnames(data)%in%preds]
preds <- colnames(data)

#----------------------------#
#### 2. SUITABILITY MODEL ####
#----------------------------#

form <- paste(preds, collapse = "+")
subs_mod <- sample(1:nrow(lu), 50000)
suitmod <- suitmodel(form = form, lu = lu[subs_mod,], data = data[subs_mod,], resolution = 10000, model = FALSE, maxit = 10000, decay = 0.01)
sm <- predict(suitmod, newdata = data, type = "probs")

mod.loglik <- nnet:::logLik.multinom(suitmod)
mod0 <- suitmodel(form = 1, data = data[subs_mod,], resolution = 10000, lu = lu[subs_mod,], decay = 0.001)
mod0.loglik <- nnet:::logLik.multinom(mod0)
mod.mfr2 <- as.numeric(1 - mod.loglik/mod0.loglik)
#Determine how much we can allocate into cells that are 0
#Turn lu data into a list of matrices (need to change code so it's stored this way to begin with)

K <- ncol(lu)

ts_inds <- 1:length(ts)
lu_obs <- list()
start <- seq(1,ncol(lu_all), by = K)

for(i in 1:length(start)){
  lu_obs[[i]] <- as.data.frame(lu_all[,start[i]:(start[i]+(K-1))])
  names(lu_obs[[i]]) <- paste0("lu", 1:K, ".obs")
  print(i)
}

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

#-----------------#
#### 3. DEMAND ####
#-----------------#

demands <- demand(landuse = lu_all, ts = ts, k = K, type = "mean")[,1:(K+1)]
demands[,1+K] <- rep(demands[1,1+K], nrow(demands)) # we assume lu 12 doesn't chnage (it's water courses)
demands[,-1] <- demands[,-1]/rowSums(demands[,-1]) #rescaling this way will change lu 12 demand again, but it's minimal so negligible.
saveRDS(demands, file = file.path(data_path, "demands.rds"))

#---------------------#
#### 4. SIMULATIONS####
#---------------------#

# Simulation parameters
params <- list(
  max_dev = 1,
  resolution = 1000000,
  growth = ch_ma[,7],
  no_change = c(9)
)

# 4.a Fully Naive model simulation

#Demand
dmd <- demands
dmd <- dmd[which(dmd[,1]%in%ts),-1]

#Prepare variables
lu_out <- lu
lu_ts <- list()
lu_suit <- list()
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
                        constraint = FALSE,
                        pa = NULL)
  
  cat('\n')
  lu_ts[[i]] <- lu_out <- lu_pred
  lu_suit[[i]] <- sm
}

#check that all predictions are within parameters
all(abs((do.call("rbind", lapply(lu_ts, FUN = colMeans)) - dmd[-1,])/dmd[-1,] * 100) < 1) #Must be TRUE
inds_pa <- which(pa == 0)
all(unlist(lapply(lu_ts, FUN = function(x) {all(x[inds_pa,]== lu[inds_pa,])}))) #must be FALSE if pa = NULL

saveRDS(lu_ts, file = file.path("outputs", "preds_naive.rds"))
saveRDS(lu_suit, file = file.path("outputs", "suit_naive.rds"))

# 4.b semi-naive model

#Demand
dmd <- demands
dmd <- dmd[which(dmd[,1]%in%ts),-1]

#Prepare variables
lu_out <- lu
lu_ts <- list()
lu_suit <- list()
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
                        constraint = FALSE,
                        pa = NULL)
  
  cat('\n')
  lu_ts[[i]] <- lu_out <- lu_pred
  lu_suit[[i]] <- sm
}

#check that all predictions are within parameters
all(abs((do.call("rbind", lapply(lu_ts, FUN = colMeans)) - dmd[-1,])/dmd[-1,] * 100) < 1) #Must be TRUE
inds_pa <- which(pa == 0)
all(unlist(lapply(lu_ts, FUN = function(x) {all(x[inds_pa,]== lu[inds_pa,])}))) #must be FALSE if pa = NULL

saveRDS(lu_ts, file = file.path("outputs", "preds_semi.rds"))
saveRDS(lu_suit, file = file.path("outputs", "suit_semi.rds"))

#4. c Full model

#Demand
dmd <- demands
dmd <- dmd[which(dmd[,1]%in%ts),-1]

#Prepare variables
lu_out <- lu
lu_ts <- list()
lu_suit <- list()
lu_pred <- matrix(NA, nrow(lu), ncol(lu))
colnames(lu_pred) <- colnames(lu)
#debug(allocation)
#Simulate time series
for(i in 1:(length(ts)-1)){
  cat(paste0('\n', "Predicting suitability for time step ", i))
  
  #Predict suitability model
  ln <- neighbourhood(lu_out, 1:K, weights, mask, enr = TRUE)
  sm <- predict(suitmod, newdata = cbind(dat, scale(ln, center = cent, scale = scal)), type = "probs")
  # sm <- matrix(runif(K*nrow(lu)), ncol = K, nrow = nrow(lu))
  # sm <- sm/rowSums(sm)
  
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
                        constraint = TRUE,
                        pa = pa)
  
  cat('\n')
  lu_ts[[i]] <- lu_out <- lu_pred
  lu_suit[[i]] <- sm
}

all(abs((do.call("rbind", lapply(lu_ts, FUN = colMeans)) - dmd[-1,])/dmd[-1,] * 100) < 1) #Must be TRUE
inds_pa <- which(pa == 0)
all(unlist(lapply(lu_ts, FUN = function(x) {all(x[inds_pa,]== lu[inds_pa,])}))) #must be TRUE if pa = NULL

saveRDS(lu_ts, file = file.path("outputs", "preds_full.rds"))
saveRDS(lu_suit, file = file.path("outputs", "suit_full.rds"))





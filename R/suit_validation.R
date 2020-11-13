# Crossvalid


rm(list = ls())

require(raster)
require(extraDistr)
require(wrswoR)
require(flutes)

source("./R/functions.R")

#---------------------------#
#### 1. DATA PREPERATION ####
#---------------------------#

# 1.a Load data
data_paths <- c(file.path(getwd(), "data", "data_ama"), file.path(getwd(), "data", "data_ama_1k"))
validation_results <- list()

for(z in 1:2){
  
  data_path <- data_paths[[z]]
  
  dat <- readRDS(file.path(data_path, "cov.rds")) #dynamic bioclimatic variables head(dat)
  lu_all <- readRDS(file.path(data_path, "lu.rds"))
  mask <- readRDS(file.path(data_path, "mask_ama.rds")) #country mask
  inds <- which(!is.na(mask[]))
  
  ts <- 1991 + c(1, 5, 10, 15, 20, 25, 27)
  K <- ncol(lu_all)/length(ts)
  n <- nrow(lu_all)
  lu <- lu_all[,1:K]
  rm(lu_all)
  
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
  
  corre <- correlations(data, sub = nrow(data)/50)
  preds <- colnames(corre)
  data <- data[,colnames(data)%in%preds]
  preds <- colnames(data)
  
  
  #----------------------------#
  #### 2. SUITABILITY MODEL ####
  #----------------------------#
  tras <- raster(crs = crs(mask), ext = extent(mask), nrow = 3, ncol = 3)
  tras[] <- 1:length(tras)
  tras <- projectRaster(tras, mask, method = "ngb")
  tras[is.na(mask)] <- NA
  
  inds <- which(!is.na(mask[]))
  tras <- tras[inds]
  
  folds <- sort(na.omit(unique(tras)))
  
  #Suitability model cross valid
  tsub <- sample(1:nrow(data), round(nrow(data)/50))
  
  data_subs <- data[tsub,]
  lu_subs <- lu[tsub,]
  tras_subs <- tras[tsub]
  
  # Build model
  partype <- c("env", "neigh", "both")
  
  # Set up model forumals for differnet models
  forms <- list(paste(preds[-which(grepl(preds, pattern = "ts"))], collapse = "+"),
                paste(preds[which(grepl(preds, pattern = "ts"))], collapse = "+"),
                paste(preds, collapse = "+"))
  
  mod <- list()
  rmse <- list()
  
  # Loop through covariate sets ("models")
  for (j in 1:length(partype)){
    rmse.mod <- data.frame("rmse" = 1:length(inds), "fold" = NA)
    coefs.mod <- list()
    # Loop through folds
    for(i in folds){
      test <- i
      test_inds <- which(tras%in%test)
      train <- folds[!folds%in%test]
      train_inds <- which(tras_subs%in%train)
      suitmod <- suitmodel(form = forms[[j]], lu = lu_subs[train_inds,], data = data_subs[train_inds,], resolution = 1000, model = FALSE, maxit = 1000, decay = 0)
      pred <- predict(suitmod, newdata = data[test_inds,], type = "probs")
      rmse.mod[test_inds, 1] <- sqrt(rowMeans((pred - lu[test_inds,])^2))
      rmse.mod[test_inds, 2] <- i
      coefs.mod[[i]] <- coefficients(suitmod)
      print(i)
    }
    rmse[[j]] <- rmse.mod
    mod[[j]] <- coefs.mod
  }
  validation_results[[z]] <- list(rmse, mod, corre)
}

saveRDS(validation_results, "validation_results.rds")

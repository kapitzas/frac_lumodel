#rm(list = ls())

require(raster)
require(matrixStats)
require(tidyverse)
require(rgdal)


source("./R/functions.R")
data_path <- file.path(getwd(), "data", "data_ama")
lu_all <- readRDS(file.path(data_path, "lu.rds"))
mask <- readRDS(file.path(data_path, "mask_ama.rds")) #country mask

preds_full <- readRDS(file.path("outputs", "preds_full.rds"))
preds_semi <- readRDS(file.path("outputs", "preds_semi.rds"))
preds_naive <- readRDS(file.path("outputs", "preds_naive.rds"))

suit_full <- readRDS(file.path("outputs", "suit_full.rds"))
suit_semi <- readRDS(file.path("outputs", "suit_semi.rds"))
suit_naive <- readRDS(file.path("outputs", "suit_naive.rds"))

mask <- readRDS(file.path(data_path, "mask_ama.rds")) #country mask

lu_all <- readRDS(file.path(data_path, "lu.rds"))

#Calculate RMSE
K <- 11
ts_inds <- 1:length(ts)
lu_obs <- list()
start <- seq(1,ncol(lu_all), by = K)
i <- 1
for(i in 1:length(start)){
  lu_obs[[i]] <- lu_all[,start[i]:(start[i]+(K-1))]
  colnames(lu_obs[[i]]) <- paste0("lu", 1:K, ".obs")
  print(i)
}

q1 <- lu_obs[[1]]
q <- lu_obs[-1]

#difference between first and consequtive time steps
q.diff <- lapply(q, FUN = function(x){rowMaxs(as.matrix(abs(x-q1)))})

#1. Null model
p_sn <- q1
null_rmse <- lapply(X = q, FUN = function(x) sqrt(rowMeans((p_sn - x)^2)))
null_mae <- lapply(X = q, FUN = function(x) rowSums(abs(p_sn - x)))
null_rsq <- lapply(X = q, FUN = function(x) {1-(rowSums((p_sn - x)^2)/(rowSums((x - colMeans(x))^2)))})

#2. RMSE
full_rmse <- list()
semi_rmse <- list()
naive_rmse <- list()
fulls_rmse <- list()
semis_rmse <- list()
naives_rmse <- list()
for(i in 1:length(q)){
  full_rmse[[i]] <- sqrt(rowMeans((preds_full[[i]] - q[[i]])^2))
  semi_rmse[[i]] <- sqrt(rowMeans((preds_semi[[i]] - q[[i]])^2))
  naive_rmse[[i]] <- sqrt(rowMeans((preds_naive[[i]] - q[[i]])^2))
  fulls_rmse[[i]] <- sqrt(rowMeans((suit_full[[i]] - q[[i]])^2))
  semis_rmse[[i]] <- sqrt(rowMeans((suit_semi[[i]] - q[[i]])^2))
  naives_rmse[[i]] <- sqrt(rowMeans((suit_naive[[i]] - q[[i]])^2))
  print(i)
}

#MAE
full_mae <- list()
semi_mae <- list()
naive_mae <- list()

for(i in 1:length(q)){
  full_mae[[i]] <- sqrt(rowSums((preds_full[[i]] - q[[i]])^2))
  semi_mae[[i]] <- sqrt(rowSums((preds_semi[[i]] - q[[i]])^2))
  naive_mae[[i]] <- sqrt(rowSums((preds_naive[[i]] - q[[i]])^2))
  print(i)
}

# full_rsq <- list()
# for(i in 1:length(q)){
#   full_rsq[[i]] <- 1-(rowSums((preds_full[[i]] - q[[i]])^2)/(rowSums((q[[i]] - colMeans(q[[i]]))^2)))
#   print(i)
#ts <- 1991 + c(5, 10, 15, 20, 25, 27)
ts <- c(5, 10, 15, 20, 25, 27)

#Make data frames
data_f1 <- cbind(tibble(
  "null_rmse" = null_rmse %>% unlist(), 
  "full_rmse" = full_rmse %>% unlist(), 
  "semi_rmse" = semi_rmse %>% unlist(), 
  "naive_rmse" = naive_rmse %>% unlist(), 
  "fulls_rmse" = full_rmse %>% unlist(), 
  "semis_rmse" = semi_rmse %>% unlist(), 
  "naives_rmse" = naive_rmse %>% unlist(), 
  "suit + constraints" = full_rmse - null_rmse,
  "suit" = semi_rmse - null_rmse,
  "naive" = naive_rmse - null_rmse,
  "year" = rep(ts, each = nrow(q1))),
  "max.diff" = q.diff %>% unlist()) %>%
  #filter(null_mae != 0) %>%
  #mutate("category" = cut(null_rmse, breaks=3))
  mutate("category" = cut(max.diff, breaks=c(-Inf, 0.005, 0.3, 0.5, Inf), 
                          labels=c("[0, 0.005)", "[0.005, 0.3)", "[0.3, 0.5)", "[0.5, inf)")))
  # mutate("category" = cut(max.diff, breaks=c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, Inf), 
  #                       labels=c("<0.1","0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "> 0.7")))
#mutate("category" = cut(null_rmse, breaks=c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, Inf), 
 #                       labels=c("<0.1","0.2", "0.3", "0.4","0.5",  ">0.5")))

#data_f2 <- data_f1 %>% gather(key = "model", value = "diffs", "diff_fulls_rmse", "diff_semis_rmse", "diff_naives_rmse")
data_f1 <- data_f1 %>% gather(key = "model", value = "diffs", "suit + constraints", "suit", "naive")
data_f1 <- data_f1 %>% gather(key = "model", value = "diffs", "null_rmse", "full_rmse", "semi_rmse", "naive_rmse")

data_aggr <- summarySE(data_f1, "diffs", groupvars = c("year", "category", "model"), conf.interval=0.95)

# Use a consistent y range
gpl <- ggplot(data_aggr, aes(x=year, y=diffs, colour=model)) + scale_colour_viridis_d(option = "viridis", begin = 0.1, end = 0.9)
gpl + geom_point() + geom_errorbar(aes(ymin=diffs-ci, ymax=diffs+ci), width=.4) + geom_line() + facet_wrap(~category) + geom_hline(yintercept = 0, linetype="dashed") + theme_bw() + scale_x_continuous(breaks = seq(5, 25, by = 5)) + ylab("RMSE difference to null")


sp <- ggplot(data_f1, aes(x = timestep, y =  diffs, fill = model))
sp + geom_boxplot(outlier.shape=".") + facet_wrap(~category) + labs(x = "Max diff", y = "RMSE difference (models - null")


#AUC DO KAPPA INSTEAD
require(pROC)
obs <- lu_obs
preds <- preds_null
aucs <- function(obs, preds){
  preds <- c(obs[1], preds[1:6])
  auc_store <- matrix(NA, ncol = 11, nrow = 6)
  roc_store  <- list()
  for (i in 2:(length(obs))){
    pred <- (preds[[i-1]] != preds[[i]]) * 1
    ob <- (obs[[i-1]] != obs[[i]]) * 1
    roc_store[[i-1]] <- list()
    j <- 1
    for (j in 1:ncol(pred)){
      roc_obj <- roc(as.numeric(ob[, j]), as.numeric(pred[, j]))
      auc_store[i-1, j] <- as.numeric(auc(roc_obj))
      roc_store[[i-1]][[j]] <- roc_obj
    }
  }
  colnames(auc_store) <- paste0("lu", str_pad(1:11, 2, pad = "0"))
  list(data.frame(auc_store), roc_store)
}

#When observed change is always 0, the true negative is 1 and true positive 0
preds_null <- list()
for(i in 1:length(preds_full)){
  preds_null[[i]] <- p_sn
}

aucs_null <- aucs(lu_obs, preds_null)
aucs_naive <- aucs(lu_obs, preds_naive)
aucs_semi <- aucs(lu_obs, preds_semi)
aucs_full <- aucs(lu_obs, preds_full)

auc_df <- rbind(aucs_naive[[1]], aucs_semi[[1]], aucs_full[[1]])
auc_df$ts <- rep(paste0("ts", 1:6), 3)
auc_df$model <- c(rep("naive", 6), rep("semi", 6), rep("full", 6))
auc_df <- gather(auc_df, key = "landuse", value = "auc", paste0("lu", str_pad(1:11, 2, pad = "0")))

sp <- ggplot(auc_df, aes(x = landuse, y =  auc, fill = model, colour = model))
sp + geom_point(position = position_dodge(0.5), shape = 19, size = 1) + theme_bw()

sp <- ggplot(auc_df, aes(x = ts, y =  auc, fill = model, colour = model))
sp + geom_point(position = position_dodge(0.5), shape = 19, size = 1) + theme_bw()

#Suitability on changing cells.
pred_store <- ob_store <- list()
obs <- lu_obs
det_changes <- function(obs){
  pred <- ob <- list()
  i <- 2
  for (i in 2:(length(obs))){
    ob[[i-1]] <- sign(obs[[i]] - obs[[i-1]])
  }
  ob
}

ch <- det_changes(c(lu_obs[1], preds_naive))

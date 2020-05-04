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
for(i in 1:length(q)){
  full_rmse[[i]] <- sqrt(rowMeans((preds_full[[i]] - q[[i]])^2))
  semi_rmse[[i]] <- sqrt(rowMeans((preds_semi[[i]] - q[[i]])^2))
  naive_rmse[[i]] <- sqrt(rowMeans((preds_naive[[i]] - q[[i]])^2))
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

#Make data frames
data_f1 <- cbind(tibble(
  "null_rmse" = null_rmse %>% unlist(), 
  "full_rmse" = full_rmse %>% unlist(), 
  "semi_rmse" = semi_rmse %>% unlist(), 
  "naive_rmse" = naive_rmse %>% unlist(), 
  "null_mae" = null_mae %>% unlist(), 
  "full_mae" = full_mae %>% unlist(), 
  "semi_mae" = semi_mae %>% unlist(), 
  "naive_mae" = naive_mae %>% unlist(), 
  "diff_full_rmse" = full_rmse - null_rmse,
  "diff_semi_rmse" = semi_rmse - null_rmse,
  "diff_naive_rmse" = naive_rmse - null_rmse,
  "Suit Model + Constraints" = full_mae - null_mae,
  "Suit Model" = semi_mae - null_mae,
  "Naive" = naive_mae - null_mae,
  "timestep" = paste0("timestep ", rep(1:length(q), each = nrow(q1)))),
  "max.diff" = q.diff %>% unlist()) %>%
  #filter(null_mae != 0) %>%
  mutate("category" = cut(max.diff, breaks=6))

#mutate("category" = cut(null_mae, breaks=c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, Inf), 
#                        labels=c("0.1","0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "> 0.7")))

data_f1 <- data_f1 %>% gather(key = "model", value = "diffs", "Suit Model + Constraints", "Suit Model", "Naive")

sp <- ggplot(data_f1, aes(x = category, y =  diffs, fill = model))
sp + geom_boxplot(position = position_dodge(0.8), outlier.shape=".") + facet_wrap(~timestep) + labs(x = "Max diff", y = "MAE difference (models - null")

sp + geom_boxplot(position = position_dodge(0.8), outlier.shape=NA) + facet_wrap(~timestep) + labs(x = "Max diff", y = "MAE difference (models - null") + ylim(-1, 0.1)

#AUC
aucs <- function(obs, preds){
  preds <- c(obs[1], preds[1:6])
  auc_store <- matrix(NA, ncol = 11, nrow = 6)
  for (i in 2:(length(obs))){
    pred <- (preds[[i-1]] != preds[[i]]) * 1
    ob <- (obs[[i-1]] != obs[[i]]) * 1
    
    for (j in 1:ncol(pred)){
      auc_store[i-1, j] <- tryCatch({
      roc_obj <- roc(as.numeric(ob[, j]), as.numeric(pred[, j]))
      as.numeric(auc(roc_obj))
      }, error=function(e) NA)
    }
  }
  colnames(auc_store) <- paste0("lu", str_pad(1:11, 2, pad = "0"))
  data.frame(auc_store)
}

#When observed change is always 0, the true negative is 1 and true positive 0
preds_null <- list()
for(i in 1:length(preds_full)){
  preds_null[[i]] <- p_sn
}

for(j in 1:length(lu_obs)){
  for(i in 1:13){
    print(length(which(lu_obs[[j]][,i] > 0)))
    #plot(makeRaster(mask, lu_obs, ts = 1, class = i))
  }
}
dev.off()


aucs_null <- aucs(lu_obs, preds_null)
aucs_naive <- aucs(lu_obs, preds_naive)
aucs_semi <- aucs(lu_obs, preds_semi)
aucs_full <- aucs(lu_obs, preds_full)

auc_df <- rbind(aucs_naive, aucs_semi, aucs_full)
auc_df$ts <- rep(paste0("ts", 1:6), 3)
auc_df$model <- c(rep("naive", 6), rep("semi", 6), rep("full", 6))
auc_df <- gather(auc_df, key = "landuse", value = "auc", paste0("lu", str_pad(1:11, 2, pad = "0")))

sp <- ggplot(auc_df, aes(x = model, y =  auc, colour = ts))
sp + geom_point() + facet_wrap(~landuse) + scale_color_brewer(palette= "Blues")
?scale_color_brewer

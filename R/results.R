rm(list = ls())

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

lu_all
#Calculate RMSE
K <- 12
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
  "full" = full_rmse - null_rmse,
  "semi-naive" = semi_rmse - null_rmse,
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


data_f1 <- data_f1 %>% gather(key = "model", value = "diffs", "full", "semi-naive", "naive")

#data_f1 <- data_f1 %>% gather(key = "model", value = "diffs", "null_rmse", "full_rmse", "semi_rmse", "naive_rmse")
colnames(data_f1)
debug(summarySE)
data_aggr <- summarySE(data_f1, measurevar = "diffs", groupvars = c("year", "category", "model"), conf.interval=0.95)

# Use a consistent y range
figure_path <- "/Users/simon/OneDrive - The University of Melbourne/PhD/writing/papers/MEE/figures/"
gpl <- ggplot(data_aggr, aes(x=year, y=diffs, colour=model)) + scale_colour_viridis_d(option = "viridis", begin = 0.1, end = 0.9)
gpl + geom_point() + 
  geom_errorbar(aes(ymin=diffs-ci, ymax=diffs+ci), width=.4) + 
  geom_line() + 
  facet_wrap(~category) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  theme_bw() + 
  scale_x_continuous(breaks = seq(5, 25, by = 5)) + 
  ylab("RMSE difference to null")

ggsave("figure2.pdf", path = figure_path, width = 18, height = 10, unit = "cm")

# #When observed change is always 0, the true negative is 1 and true positive 0
preds_null <- list()
for(i in 1:length(preds_full)){
  preds_null[[i]] <- p_sn
}

ref <- 1
aucs_null <- diff_metrics(lu_obs, preds_null, mask, percent = TRUE, reference = ref)
aucs_naive <- diff_metrics(lu_obs, preds_naive, mask, percent = TRUE, reference = ref)
aucs_semi <- diff_metrics(lu_obs, preds_semi, mask, percent = TRUE, reference = ref)
aucs_full <- diff_metrics(lu_obs, preds_full, mask, percent = TRUE, reference = ref)
aucs <- list(aucs_null, aucs_naive, aucs_semi, aucs_full)

require(diffeR)

t1 <- do.call("rbind", lapply(aucs, FUN = function(x) 
  do.call("rbind", lapply(x, FUN = function(x) 
    do.call("c", lapply(x, FUN = function(x) 
      overallDiff(x)))))))
t1 <- as.data.frame(t1)

t1$lu <- as.factor(rep(1:K, 4))
t1$model <- factor(rep(c("null", "naive", "semi", "full"), each = 12), levels = c("null", "naive", "semi", "full"))
t1 <- t1 %>% gather(key = "timestep", value = "diff", -lu, - model)
t1$timestep <- as.factor(rep(c(5, 10, 15, 20, 25, 27), each = 48))

sp <- ggplot(t1, aes(x = timestep, y = diff, colour = model))

sp + geom_point(position=position_dodge(0.3), shape = 16) + scale_colour_viridis_d(option = "viridis", begin = 0.95, end = 0.05) + theme_bw() + ylab("disagreement [% of observations]") + xlab("Year")

t1$timestep <- as.factor(paste0("timestep ", rep(c("05", 10, 15, 20, 25, 27), each = 4)))

# sp <- ggplot(t1, aes(x = model, y = diff, colour = lu))
# sp + geom_point(position=position_dodge(0.3), shape = 16) + 
#   scale_colour_viridis_d(option = "viridis", begin = 0.95, end = 0.05) + 
#   theme_bw() + 
#   ylab("disagreement [% of observations]") + 
#   xlab("Model") + facet_wrap(~timestep) + 
#   theme(text=element_text(size=12))

ggsave("figure3.pdf", path = figure_path, width = 18, height = 10, unit = "cm", dpi = 600)

# t1$timestep <- as.factor(paste0("timestep ", rep(c("05", 10, 15, 20, 25, 27), each = 4)))
# sp <- ggplot(t1, aes(x = model, y = diff, colour = lu))
# sp + geom_point(position=position_dodge(0.3), shape = 16) + scale_colour_viridis_d(option = "viridis", begin = 0.95, end = 0.05) + theme_bw() + ylab("disagreement [% of observations]") + xlab("Model") + facet_wrap(~timestep)

#+ geom_line(aes(group = interaction(lu, model)), position = position_dodge(0.6), size = 0.1)
#Road map for the next part:

#1. Determine which land use classes are important for maybe 10 threatened species (observations) (using predicts)
predicts <- readRDS("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter2/PREDICTS/database.rds")

dom_obs <- lapply(lu_obs, FUN = function(x) { as.numeric(apply(x, 1, FUN =  function(x) { which(x == max(x), arr.ind = TRUE)[1]}))})
dom_full <- lapply(preds_full, FUN = function(x) {as.numeric(apply(x, 1, FUN =  function(x) { which(x == max(x), arr.ind = TRUE)[1]}))}) #not the same as above

diff_obs <- lu_obs[[7]] - lu_obs[[1]]
diff_full <- preds_full[[6]] - lu_obs[[1]]

chcrop_obs <- rowSums(diff_obs[,c(1,2)]) #cropland expansion per cell
#chcrop_obs <- diff_obs[,c(1)] #cropland expansion per cell
chnat_obs <- rowSums(diff_obs[,c(3,4,5, 8, 9)])
#chnat_obs <- diff_obs[,c(3)]

chcrop_full <- rowSums(diff_full[,c(1,2)]) #cropland expansion per cell
#chcrop_full <- diff_full[,c(1)] #cropland expansion per cell
chnat_full <- rowSums(diff_full[,c(3,4,5,6,8,9)])
#chnat_full <- diff_full[,c(3)]

inds <- which(!is.na(mask[]))
r <- mask
r[inds] <- 0

hloss_full <- hloss_obs <- r
hloss_full[inds[which(chcrop_full > 0 & chnat_full < 0)]] <- 1
hloss_obs[inds[which(chcrop_obs > 0 & chnat_obs < 0)]] <- 1
hloss_full[inds[which(chcrop_full < 0 & chnat_full > 0)]] <- 2
hloss_obs[inds[which(chcrop_obs < 0 & chnat_obs > 0)]] <- 2


plot(hloss_full)
plot(hloss_obs)
r1 <- r2 <- r
par(mfrow = c(7, 7))
for(i in 1: ncol(data)){
  plot(diff_obs[,2] ~ data[,i])
}
r1[inds] <- diff_obs[,2]

r2[inds] <- diff_full[,2]

r1[inds] <- lu_obs[[7]][,2]
r2[inds] <- preds_full[[6]][,2]

r3 <- makeRaster(mask, lu = sm, class = 2)
rasterVis::levelplot(stack(r1, r2))

preds
plot(makeRaster(mask, data, class = 13)
rasterVis::levelplot(stack(r1,r3))
makeRaster(sm)
plot(r)
hist(diff_full[,1])
hist(diff_obs[,1])
overallDiff(crosstabm(hloss_obs, hloss_full, percent = TRUE))
hloss_full <- hloss_obs <- r
hloss_full[inds[which(chnat_full < 0)]] <- 1
hloss_obs[inds[which(chnat_obs < 0)]] <- 1
overallDiff(ctmatrix = )
hloss_full <- hloss_obs <- r
hloss_full[inds[which(chcrop_full > 0)]] <- 1
hloss_obs[inds[which(chcrop_obs > 0)]] <- 1

plot(hloss_obs)
plot(hloss_full)
plot(makeRaster(mask, lu_obs, class = 2, ts = 1))
plot(makeRaster(mask, sm, class = 2))
store <- list()
n <- nrow(lu_obs[[1]])
thresh <- seq(0, 1, by = 0.1)
tables_crop <- tables_nat <- store_nat <- store_crop <- list()
for(i in seq_along(thresh)[-11]){
  
  hloss_full <- hloss_obs <- r
  hloss_full[inds[which(chcrop_full > thresh[i] & chcrop_full < thresh[i+1])]] <- 1
  hloss_obs[inds[which(chcrop_obs > thresh[i] & chcrop_obs < thresh[i+1])]] <- 1
  tables_crop[[i]] <- crosstabm(hloss_obs, hloss_full, percent = TRUE)
  store_crop[[i]] <- overallDiff(tables_crop[[i]])
  
  hloss_full <- hloss_obs <- r
  hloss_full[inds[which(chnat_full < -thresh[i] & chnat_full > -thresh[i+1])]] <- 1
  hloss_obs[inds[which(chnat_obs < -thresh[i] & chnat_obs > -thresh[i+1])]] <- 1
  tables_nat[[i]] <- crosstabm(hloss_obs, hloss_full, percent = TRUE)
  store_nat[[i]] <- overallDiff(tables_nat[[i]])
}

plot(thresh[-11], unlist(store_crop), type = "l", col = "blue", lwd = 1.2)
lines(thresh[-11], unlist(store_nat), col = 'red', lwd = 1.2)

r[r%in%crop] <- 1
r[r%in%crop_natveg_mosaic] <- 2
r[r%in%tree_co] <- 3
r[r%in%treeshrub_herb_mosaic] <- 4
r[r%in%shrub] <- 5
r[r%in%grass] <- 6
r[r%in%sparse] <- 7
r[r%in%tree_water] <- 8
r[r%in%shrub_water] <- 9
r[r%in%urban] <- 10
r[r%in%bare] <- 11
r[r%in%snow_ice] <- 12
r[r%in%water] <- 13



diffeR::agreementj()
which(nat_red)

plot(forest_red ~ cropland_exp)

cropland_exp_preds <- rowSums(diff_full[,c(1,2)]) #cropland expansion per cell

forest_red_preds <- rowSums(diff_full[,c(4,5)])

plot(forest_red_preds ~ cropland_exp_preds)



diff_obs[[which(diff_obs)]]


d1_obs[inds] <- dom_obs[[1]]
d2_full[inds] <- dom_full[[6]]
d2_obs[inds] <- dom_obs[[7]]

overallDiffCatj(crosstabm(d1_obs, d2_obs))
overallDiffCatj(crosstabm(d1_obs, d2_full))

overallDiffCatj(d1_obs, d2_obs)



diffeR::exchangeDij(crosstabm(d1, d2))


n <- nrow(lu_obs[[1]])



#2. Aggregate our our mapping to those classes (observations + simulatons)
#3. Classify to dominant land use (observations + simulations)
#4. Determine transitions from natural classes to cropland class/classes (observations + simulations)
#5. classify all maps to dominant land use classes
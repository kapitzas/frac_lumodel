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
K <- 13
ts_inds <- 1:length(ts)
lu_obs <- list()
start <- seq(1,ncol(lu_all), by = K)
for(i in 1:length(start)){
  lu_obs[[i]] <- as.data.frame(lu_all[,start[i]:(start[i]+(K-1))])
  names(lu_obs[[i]]) <- paste0("lu", 1:K, ".obs")
  print(i)
}

q1 <- lu_obs[[1]]
q <- lu_obs[-1]

rm(lu_all, lu_obs)

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
  "diff_full_mae" = full_mae - null_mae,
  "diff_semi_mae" = semi_mae - null_mae,
  "diff_naive_mae" = naive_mae - null_mae,
  "timestep" = rep(1:length(q), each = nrow(q1))),
  "max.diff" = q.diff %>% unlist()) %>%
  mutate("category" = cut(null_mae, breaks=c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, Inf), 
                          labels=c("0.1","0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "> 0.9")))

data_f1 <- data_f1 %>% gather(key = "model", value = "diffs", diff_full_mae, diff_semi_mae, diff_naive_mae)

sp <- ggplot(data_f1, 
             aes(x = category, y =  diffs, fill = model))


sp + geom_boxplot(position = position_dodge(0.8), outlier.shape=".") + facet_wrap(~timestep) + labs(x = "Levels of change", y = "Difference MAE")


sp + geom_boxplot(position = position_dodge(0.8), outlier.shape=NA) + facet_wrap(~timestep) + labs(x = "Levels of change", y = "Difference RMSE") + ylim(-1, 0.1)


sp <- ggplot(data_f2, 
             aes(x = category, y =  diffs, fill = model))

sp <- sp + geom_boxplot(outlier.shape = ".")  + facet_wrap(~timestep) + ylim(-1, 1)
sp <- sp + geom_boxplot(outlier.shape = ".")  + facet_wrap(~timestep)

sp + labs(y="RMSE", x = "Maximum change in any class from first time step")


preds_cat <-  preds_full %>% map(as.data.frame) %>% bind_rows()
plot_data <- boxplot(formula = diffs ~ model + category, data = data_f1, outline = FALSE, plot = FALSE)
plot_data <- boxplot(formula = diffs ~ model + category, data = data_f1, outline = FALSE, plot = TRUE)
no_outliersn <- which(!1:max(unique(plot_data$group))%in%unique(plot_data$group))
not_plotted <- round(table(out$group)/out$n[-no_outliers] * 100, 2)
not_plotted <- c(not_plotted[c(1:18)], "19" = 0,  not_plotted[c(19)])



make_breaks <- function(x) {
  cut(x, breaks=c(-Inf, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, Inf), 
      labels=c("0.05","0.1","0.2", "0.3", "0.4", "0.5", "0.5"))
}

preds_cat <- as.data.frame(lapply(preds_cat, FUN = make_breaks))

data_f1 <- cbind(data_f1, preds_cat)
starts_with()
data_f2 <- data_f1 %>% sample_n(1000000)
data_f2 <- data_f2 %>% gather(key = "model", value = "diffs", diff_full, diff_random)
data_f2 <- data_f2 %>% gather(key = "class", value = "level", starts_with("ts"))
str(data_f2)
mask

cl <- 2
par(mfrow = c(1, 2))
plot(makeRaster(mask, lu = preds_full, ts = 6, class = cl))
plot(makeRaster(mask, lu = q, ts = 1, class = cl))

(preds_full[[6]][,11])
sp <- ggplot(data_f2 %>% sample_n(5000000), 
             aes(x = level, y =  diffs, fill = model))
sp <- sp + geom_boxplot(position = position_dodge(0.8), outlier.size = 0)
sp + facet_wrap(~class)

sp + labs(y="RMSE", x = "Maximum change in any class from first timestep")





sp <- ggplot(data_f1 %>% sample_n(500000), 
             aes(x = category, y =  diffs, fill = model))

sp <- sp + geom_boxplot(position = position_dodge(0.8), outlier.size = 0)
sp + labs(y="RMSE", x = "Maximum change in any class from first timestep")


# No changes pl

sp <- ggplot(data_f1 %>% sample_n(size = 500000),
             aes(x = max.diff, y =  diff_rmse))
sp + geom_point()
sp <- sp + geom_boxplot(position = position_dodge(0.8), outlier.size = 0)
sp



n <- nrow(lu_ts[[1]])

# Plots

#1
sp <- ggplot(data_f1 %>% sample_n(size = 500000),
             aes(x = timestep, y =  log(rmse), fill = model))
sp <- sp + geom_boxplot(position = position_dodge(0.8), outlier.size = 0)
sp

#sp <- sp + facet_wrap(~country)

ggsave("figure2.pdf", path = "/Users/simon/OneDrive - The University of Melbourne/PhD/writing/papers/MEE/figures/")

#2
# sp <- ggplot(data_f2 %>% sample_n(size = 500000),
#              aes(x = timestep, y = difference, fill = class))
# sp <- sp + geom_boxplot(position = position_dodge(0.8), outlier.size=0)
# sp <- sp + facet_wrap(~country)
# ggsave("figuren.pdf", path = "/Users/simon/OneDrive - The University of Melbourne/PhD/writing/papers/MEE/figures/")
#3

## Use densCols() output to get density at each point
# countries <- unique(data_f3$country)
# cols <- viridis::plasma(256)
# dens <- numeric()
# col <- numeric()
# for(i in 1:length(countries)){
#   sub <- data_f3 %>% 
#     filter(country == countries[i]) %>%
#     select(null, diff_rmse)
#   x <- densCols(sub$null, sub$diff_rmse, colramp=colorRampPalette(c("black", "white")))
#   dens <- c(dens, col2rgb(x)[1,] + 1L)
#   col <- c(dens, cols[dens])
#   print(i)
# }
# plot(diff_rmse ~ null, data=td[order(td$dens),], pch= 17, col=col)

## Plot it, reordering rows so that densest points are plotted on top


sp <- sp + stat_density2d(aes(fill = ..density..^0.1), geom = "tile", contour = FALSE, n = 2000) +   scale_fill_continuous(low = "white", high = "dodgerblue4") + geom_hline(yintercept=0, linetype="solid", color = "darkgreen", size=0.5)
#qplot(test_data$null, diff_rmse, data = augment(fit), col = alpha(0.5)) + geom_line(aes(y = .fitted), color = "darkgreen")

sp <- ggplot(data_f3 %>% sample_n(size = 100000),
             aes(x = null, y = diff_rmse))
sp <- sp + geom_point(position=position_jitter(h=0.01, w=0.01), size = 0.01, alpha = 0.3)
sp <- sp + geom_hline(yintercept=0, linetype="solid", 
                      color = "red", size=0.5)
sp +  geom_smooth(method = "lm", formula = y ~ poly(x, 3), size = 0.5)
sp


ggsave("figure3.pdf", path = "/Users/simon/OneDrive - The University of Melbourne/PhD/writing/papers/MEE/figures/")

#4
# sp <- ggplot(data_f2 %>% sample_n(size = 20000),
#              aes(x = difference, y = rmse, color = timestep))
# sp <- sp + geom_point()
# sp + facet_wrap(~class+timestep)

test <- lu_obs[[2]] - lu_ts[[1]]
cor(as.numeric(lu_ts[[1]][1,]), as.numeric(lu_obs[[2]][1,]))
cor_full <- cor_null <- numeric()
for(i in 1:nrow(lu_obs[[1]])){
  cor_null[i] <- cor(as.numeric(lu_obs[[1]][i,]), as.numeric(lu_obs[[2]][i,]))
  cor_full[i] <- cor(as.numeric(lu_ts[[1]][i,]), as.numeric(lu_obs[[2]][i,]))
  print(i)
}

plot(cor_full-cor_null ~cor_null)
abline(0,0, 1, 1)
plot(makeRaster(mask, test, class = 5))
?cor
round(diag(cor(lu_ts[[1]], lu_obs[[2]])),3)
round(diag(cor(lu_ts[[2]], lu_obs[[3]])),3)
round(diag(cor(lu_ts[[3]], lu_obs[[4]])),3)

round(diag(cor(lu_obs[[1]], lu_obs[[2]])),3)
round(diag(cor(lu_obs[[1]], lu_obs[[4]])),3)
?cor

which(is.infinite(p[,4]))

rmse_map <- stack(mask, mask, mask)
for(i in 1:3){
  rmse_map[[i]][!is.na(mask[])] <- full_whole[[i]]
  print(i)
}
names(rmse_map) <- paste0("TS", 1:3)
png("/Users/simon/OneDrive - The University of Melbourne/PhD/writing/papers/MEE/figures/figure4.png")
print(levelplot(rmse_map, margin = F))
dev.off()

require(reshape2)
length(lu_obs)
obs_diff <- lu_obs[[4]] - lu_obs[[1]]
colnames(obs_diff) <- paste0("lu_", 1:5, "_diff")
observations <- melt(obs_diff)

envdata <- melt(cbind(data, ln))
envdata <- rbind(envdata[,c(2,3)], observations[,c(2,3)])
test <- data.frame("FULL_RMSE" = plot_data$full[which(plot_data$timesteps == 1)], envdata, 
                   "NULL_RMSE" = plot_data$sn[which(plot_data$timesteps == 1)], 
                   "DIFF" = plot_data$full[which(plot_data$timesteps == 1)] - plot_data$sn[which(plot_data$timesteps == 1)])
test$Var2 <- factor(test$Var2, levels = as.character(unique(test$Var2)))

test <- test[sample(1:nrow(test), 100000),]

ggplot(test, aes(y=FULL_RMSE, x=value)) + geom_point(shape=20, size = 1) +  facet_wrap(~Var2, scales = "free")

max(lu_obs[[1]][,1] - lu_ts[[3]][,1]) # 0.5424 without constraint
test <- mask
obs <- makeRaster(mask, lu = lu_obs, ts = 1, class = 1)
pred <- makeRaster(mask, lu = lu_ts, ts = 1, class = 1)

plot(obs)

#rm(list = ls())

require(raster)
require(matrixStats)
require(tidyverse)
require(rgdal)
require(rgeos)

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
  "semi" = semi_rmse - null_rmse,
  "naive" = naive_rmse - null_rmse,
  "Year" = rep(ts, each = nrow(q1))),
  "max.diff" = q.diff %>% unlist()) %>%
  #filter(null_mae != 0) %>%
  #mutate("category" = cut(null_rmse, breaks=3))
  mutate("category" = cut(max.diff, breaks=c(-Inf, 0.005, 0.3, 0.5, Inf), 
                          labels=c("[0, 0.005)", "[0.005, 0.3)", "[0.3, 0.5)", "[0.5, inf)")))


# mutate("category" = cut(max.diff, breaks=c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, Inf),
#                       labels=c("<0.1","0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "> 0.7")))
# mutate("category" = cut(null_rmse, breaks=c(-Inf, 0.1, 0.2, 0.3, 0.4, 0.5, Inf),
#                       labels=c("<0.1","0.2", "0.3", "0.4","0.5",  ">0.5")))


data_f1 <- data_f1 %>% gather(key = "model", value = "diffs", "full", "semi", "naive")

#data_f1 <- data_f1 %>% gather(key = "model", value = "diffs", "null_rmse", "full_rmse", "semi_rmse", "naive_rmse")

col.positions <- c(1,2,3,4)
require(scales)
rescale(col.positions)
data_aggr <- summarySE(data_f1, measurevar = "diffs", groupvars = c("Year", "category", "model"), conf.interval=0.95)

# Use a consistent y range
data_aggr$model <- factor(data_aggr$model, c("full", "semi", "naive"))
figure_path <- "/Users/simon/OneDrive - The University of Melbourne/PhD/writing/papers/MEE/figures/"
fig2 <- ggplot(data_aggr, aes(x=Year, y=diffs, colour=model)) + scale_colour_viridis_d(option = "plasma", begin = 0.2, end = 0.6)
fig2 <- fig2 + geom_point(shape = 16, size = 2) + 
  geom_errorbar(aes(ymin=diffs-ci, ymax=diffs+ci), width=.4) + 
  geom_line() + 
  facet_wrap(~category, scales='free') + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  theme_bw() + 
  scale_x_continuous(breaks = seq(5, 27, by = 5)) + 
  scale_y_continuous(limits=c(-0.06, 0.03)) +
  ylab("RMSE difference to null") + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white"),
        strip.text = element_text(size=10, colour="black"), 
        #legend.position = c(.5, .65), 
        legend.text = element_text(size = 10, colour = "black"),
        legend.title =  element_blank())


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

t1$model <- factor(t1$model, c("full", "semi", "naive", "null"))
fig2b <- ggplot(t1, aes(x = timestep, y = diff, colour = model))

fig2b <- fig2b + geom_point(position=position_dodge(0.5), shape = 16, size = 2) + 
  scale_colour_viridis_d(option = "plasma", begin = 0.2, end = 0.8) + 
  theme_bw() + 
  ylab("disagreement [% of observations]") + 
  xlab("Year") +
  theme(
    legend.position = "right", 
    legend.title =  element_blank(),
    legend.text = element_text(size = 10, colour = "black"),
    legend.margin=margin(c(0,0,0,0)))

plot_grid(fig2, fig2b, nrow = 2, labels = c("a", "b"), vjust = c(1.5, -1), rel_heights = c(3,2))

# sp <- ggplot(t1, aes(x = model, y = diff, colour = lu))
# sp + geom_point(position=position_dodge(0.3), shape = 16) + 
#   scale_colour_viridis_d(option = "viridis", begin = 0.95, end = 0.05) + 
#   theme_bw() + 
#   ylab("disagreement [% of observations]") + 
#   xlab("Model") + facet_wrap(~timestep) + 
#   theme(text=element_text(size=12))


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

crop1_obs <- diff_obs[,c(1)] #cropland expansion per cell
crop2_obs <- diff_obs[,c(2)] #cropland expansion per cell
nat1_obs <- diff_obs[,c(3)]
nat2_obs <- rowSums(diff_obs[,c(4,5, 8, 9)])

crop1_full <- diff_full[,c(1)] #cropland expansion per cell
crop2_full <- diff_full[,c(2)] #cropland expansion per cell
nat1_full <- diff_full[,c(3)]
nat2_full <- rowSums(diff_full[,c(4,5, 8, 9)])
hloss <- agr_
inds <- which(!is.na(mask[]))
r <- mask
r[inds] <- 0

agrexp_full <- agrexp_obs <- hloss_full <- hloss_obs <- r

agrexp_full[inds[which(rowSums(cbind(crop1_full, crop2_full)) > 0 & crop1_full > crop2_full)]] <- 1
agrexp_full[inds[which(sum(crop1_full, crop2_full) > 0 & crop2_full > crop1_full)]] <- 2

agrexp_obs[inds[which(rowSums(cbind(crop1_full, crop2_full)) > 0 & crop1_obs > crop2_obs)]] <- 1
agrexp_obs[inds[which(rowSums(cbind(crop1_full, crop2_full)) > 0 & crop2_obs > crop1_obs)]] <- 2

hloss_full[inds[which(rowSums(cbind(nat2_full, nat1_full)) < 0 & nat1_full < nat2_full)]] <- 1
hloss_full[inds[which(rowSums(cbind(nat2_full, nat1_full)) < 0 & nat2_full < nat1_full)]] <- 2

hloss_obs[inds[which(rowSums(cbind(nat2_full, nat1_full)) < 0 & nat1_obs < nat2_obs)]] <- 1
hloss_obs[inds[which(rowSums(cbind(nat2_full, nat1_full)) < 0 & nat2_obs < nat1_obs)]] <- 2

boundary <- readOGR("/Users/simon/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/Amazon boundaries/Lim_Biogeografico.shp")

?st_simplify
boundary <- gSimplify(boundary, tol = 0.05)
require(rasterVis)
require(viridis)
s <- stack(agrexp_obs, agrexp_full, hloss_obs, hloss_full)

maps <- list()
for (i in 1: nlayers(s)){
  
  r <- ratify(s[[i]])
  rat <- levels(r)[[1]]
  rat$ID <- c("0", "1", "2")
  levels(r) <- rat
  
  l <- levelplot(r, att = "ID",
                 margin=FALSE,                       
                 colorkey=FALSE,    
                 xlab="", ylab="",
                 par.settings=list(
                   strip.border=list(col='transparent'),
                   strip.background=list(col='transparent'),
                   axis.line=list(col='transparent')
                 ),
                 scales=list(draw=FALSE),            
                 col.regions=plasma(3, begin = 0.2, end = 0.8, direction = -1),
                 names.attr=rep('', nlayers(r))) +           
    layer(sp.polygons(boundary, lwd=2))
  
  if(i%in%c(1,3)){
    l <- l + layer(panel.text(-55, 9, labels = "observed", cex = 0.8, adj = 0, font = 2))
  }
  if(i%in%c(2,4)){
    l <- l + layer(panel.text(-55, 9, labels = "predicted", cex = 0.8, adj = 0, font = 2))
  }
  maps[[i]] <- l
}


legends <- list()
for(i in 1:2){
  
  if(i == 1){
    lab <- c("no increase", "dominated by cropland increase", "dominated by cropland mosaic increase")
  }
  if(i == 2){
    lab <- c("no decrease", "dominated by forest loss", "dominated by other natural habitat loss")
  }
  
  leg_data <- tibble(label = lab)
  leg_data$label <- factor(leg_data$label, lab)
  
  my_hist <- ggplot(leg_data, aes(label, fill = label)) + 
    geom_bar() + 
    scale_fill_manual(values= plasma(3, begin = 0.2, end = 0.8, direction = -1)) + 
    theme(title = element_blank(), 
          legend.position="bottom", 
          legend.margin=margin(t=-0.5, r=-0.5, b=-0.5, l=-0.5, unit="cm"))
  legends[[i]] <- cowplot::get_legend(my_hist)
}

plot_grid(arrangeGrob(maps[[1]], maps[[2]], ncol = 2), legends[[1]], arrangeGrob(maps[[3]], maps[[4]], ncol = 2), legends[[2]], nrow = 4, 
          rel_heights = c(1/1.25, 1/12, 1/1.25, 1/12), labels = c("a", "", "b", ""))


diff2 <- diff1 <- r
diff1[inds] <- rowSums(cbind(crop1_obs, crop2_obs))
diff2[inds] <- rowSums(cbind(nat1_obs, nat2_obs))
?levelplot
diffstack <- stack(stack(diff1, diff2))
p1 <- levelplot(diffstack, nrow = 1,
          margin=FALSE,                       
          colorkey=list(space='bottom'),
          xlab="", ylab="",
          par.settings=list(
            strip.border=list(col='transparent'),
            strip.background=list(col='transparent'),
            axis.line=list(col='transparent')
          ),
          scales=list(draw=FALSE),
          names.attr=rep('', nlayers(diffstack)),
          col.regions=gray(0:30/30)) +
  layer(sp.polygons(boundary, lwd=2))
p1
install.packages("stargazer")
require(viridis)
?grid.arrange
?geom_raster
require(ggspatial)
geom_raster()
levelplot(hloss_obs)
levelplot(hloss_full)
levelplot(agrexp_full)
levelplot(agrexp_obs)
plot(hloss_obs)

hloss_obs[inds[which(crop1_obs > 0 & nat1_obs < 0)]] <- 1
plot(hloss_obs)

hloss_full[inds[which(crop1_full > 0 & nat2_full < 0)]] <- 2
hloss_obs[inds[which(crop1_obs > 0 & nat2_obs < 0)]] <- 2

hloss_full[inds[which(crop2_full > 0 & nat1_full < 0)]] <- 3
hloss_obs[inds[which(crop2_obs > 0 & nat1_obs < 0)]] <- 3

hloss_full[inds[which(crop2_full > 0 & nat2_full < 0)]] <- 4
hloss_obs[inds[which(crop2_obs > 0 & nat2_obs < 0)]] <- 4


plot(hloss_full)
plot(hloss_obs)
r3 <- r4 <- r1 <- r2 <- r

r1[inds] <- diff_obs[,2]
r2[inds] <- diff_full[,2]
r3[inds] <- diff_obs[,3]
r4[inds] <- diff_full[,3]


r3 <- makeRaster(mask, lu = sm, class = 2)
rasterVis::levelplot(r2-r1)
rasterVis::levelplot(stack(r1<0, r2<0))

plot(r2)
rasterVis::levelplot(stack(r1,r2, r3, r4))
plot(r)

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
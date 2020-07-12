#rm(list = ls())

require(raster)
require(matrixStats)
require(tidyverse)
require(rgdal)
require(rgeos)
require(diffeR)
require(rasterVis)
require(viridis)
require(cowplot)
require(sf)

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
K <- 9
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




data_f1 <- data_f1 %>% gather(key = "model", value = "diffs", "full", "semi", "naive")

data_aggr <- summarySE(data_f1, measurevar = "diffs", groupvars = c("Year", "category", "model"), conf.interval=0.95)

# Use a consistent y range
data_aggr$model <- factor(data_aggr$model, c("full", "semi", "naive"))
figure_path <- "/Users/simon/OneDrive - The University of Melbourne/PhD/writing/papers/MEE/figures/"
fig3a <- ggplot(data_aggr, aes(x=Year, y=diffs, colour=model)) + scale_colour_viridis_d(option = "plasma", begin = 0.2, end = 0.6)
fig3a <- fig3a + geom_point(shape = 16, size = 2) + 
  geom_errorbar(aes(ymin=diffs-ci, ymax=diffs+ci), width=.4) + 
  geom_line() + 
  facet_wrap(~category, scales='free') + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  theme_bw() + 
  scale_x_continuous(breaks = seq(5, 27, by = 5)) + 
  scale_y_continuous(limits=c(-0.08, 0.03)) +
  ylab("RMSE difference to null") + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white"),
        strip.text = element_text(size=8, colour="black"), 
        #legend.position = c(.5, .65), 
        legend.text = element_text(size = 8, colour = "black"),
        legend.title =  element_blank())


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



t1 <- do.call("rbind", lapply(aucs, FUN = function(x) 
  do.call("rbind", lapply(x, FUN = function(x) 
    do.call("c", lapply(x, FUN = function(x) 
      overallDiff(x)))))))
t1 <- as.data.frame(t1)

t1$lu <- as.factor(rep(1:K, 4))
t1$model <- factor(rep(c("null", "naive", "semi", "full"), each = K), levels = c("null", "naive", "semi", "full"))
t1 <- t1 %>% gather(key = "timestep", value = "diff", -lu, - model)
t1$timestep <- as.factor(rep(c(5, 10, 15, 20, 25, 27), each = K * 4))

t1$model <- factor(t1$model, c("full", "semi", "naive", "null"))
fig3b <- ggplot(t1, aes(x = timestep, y = diff, colour = model))

fig3b <- fig3b + geom_point(position=position_dodge(0.5), shape = 16, size = 2) + 
  scale_colour_viridis_d(option = "plasma", begin = 0.2, end = 0.8) + 
  theme_bw() + 
  ylab("disagreement [% of cells]") + 
  xlab("Year") +
  theme(
    legend.position = "right", 
    legend.title =  element_blank(),
    legend.text = element_text(size = 8, colour = "black"),
    legend.margin=margin(c(0,0,0,0)))


png(file.path(figure_path, "figure3.png"), height = 140, width = 140, units = "mm", res = 1000)
plot_grid(fig3a, fig3b, nrow = 2, labels = c("a)", "b)"), vjust = c(1.5, -1), rel_heights = c(3,2))
dev.off()

#1. Determine which land use classes are important for maybe 10 threatened species (observations) (using predicts)
dom_obs <- lapply(lu_obs, FUN = function(x) { as.numeric(apply(x, 1, FUN =  function(x) { which(x == max(x), arr.ind = TRUE)[1]}))})
dom_full <- lapply(preds_full, FUN = function(x) {as.numeric(apply(x, 1, FUN =  function(x) { which(x == max(x), arr.ind = TRUE)[1]}))}) #not the same as above
inds <- which(!is.na(mask[]))

agrexp_hloss <- list()
i <- 1
for (i in 1:length(preds_full)){
  
  r <- mask
  r[inds] <- 0
  
  diff_obs <- lu_obs[[i+1]] - lu_obs[[1]] #observed
  diff_full <- preds_full[[i]] - lu_obs[[1]] #predicted
  
  
  crop1_obs <- diff_obs[,c(1)] #cropland expansion per cell
  crop2_obs <- diff_obs[,c(2)] #cropland expansion per cell
  pasture_obs <- diff_obs[,c(4,5)] #pasture expansion per cell
  
  nat1_obs <- diff_obs[,c(3)]
  nat2_obs <- rowSums(diff_obs[,c(6,8)])
  
  crop1_full <- diff_full[,c(1)] #cropland expansion per cell
  crop2_full <- diff_full[,c(2)] #cropland expansion per cell
  pasture_full <- diff_full[,c(4,5)] #pasture expansion per cell
  
  nat1_full <- diff_full[,c(3)]
  nat2_full <- rowSums(diff_full[,c(6,8)])
  
  
  pastureexp_obs <- pastureexp_full <- agrexp_full <- agrexp_obs <- hloss_full <- hloss_obs <- r
  
  #PASTURE EXPANSION
  #OBSERVED
  #Pasture cuts into forest
  pastureexp_obs[inds[which(pasture_obs > 0 & rowSums(cbind(nat1_obs,nat2_obs)) < 0 & nat1_obs < nat2_obs)]] <- 1
  
  #Pasture cuts into shrub/other
  pastureexp_obs[inds[which(pasture_obs > 0 & rowSums(cbind(nat1_obs,nat2_obs)) < 0 & nat2_obs < nat1_obs)]] <- 2
  
  #PREDICTED
  #Pasture cuts into forest
  pastureexp_full[inds[which(pasture_full > 0 & rowSums(cbind(nat1_full,nat2_full)) < 0 & nat1_full < nat2_full)]] <- 1
  
  #Pasture cuts into shrub/other
  pastureexp_full[inds[which(pasture_full > 0 & rowSums(cbind(nat1_full,nat2_full)) < 0 & nat2_full < nat1_full)]] <- 2
  
  #CROPLAND EXPANSION
  #PREDICTED
  #cuts into forest
  agrexp_full[inds[which(rowSums(cbind(crop1_full, crop2_full)) > 0 & rowSums(cbind(nat1_full,nat2_full)) < 0 & nat1_full < nat2_full)]] <- 1
  
  #cuts into shrub/other
  agrexp_full[inds[which(rowSums(cbind(crop1_full, crop2_full)) > 0 & rowSums(cbind(nat1_full,nat2_full)) < 0 & nat2_full < nat1_full)]] <- 2
  
  #cuts into forest
  agrexp_obs[inds[which(rowSums(cbind(crop1_obs, crop2_obs)) > 0 & rowSums(cbind(nat1_obs,nat2_obs)) < 0 & nat1_obs < nat2_obs)]] <- 1
  
  #cuts into shrub/other
  agrexp_obs[inds[which(rowSums(cbind(crop1_obs, crop2_obs)) > 0 & rowSums(cbind(nat1_obs,nat2_obs)) < 0 & nat2_obs < nat1_obs)]] <- 2
  
  # #PREDICTED
  # #Cropland dominates increase
  # agrexp_full[inds[which(rowSums(cbind(crop1_full, crop2_full)) > 0 & crop1_full > crop2_full)]] <- 1
  # 
  # #Cropland mosaic dominates increase
  # agrexp_full[inds[which(sum(crop1_full, crop2_full) > 0 & crop2_full > crop1_full)]] <- 2
  # 
  # #OBSERVED
  # #Cropland dominates increase
  # agrexp_obs[inds[which(rowSums(cbind(crop1_full, crop2_full)) > 0 & crop1_obs > crop2_obs)]] <- 1
  # 
  # #Cropland mosaic dominates increase
  # agrexp_obs[inds[which(rowSums(cbind(crop1_full, crop2_full)) > 0 & crop2_obs > crop1_obs)]] <- 2
  
  #HABITAT LOSS
  #PREDICTED
  #Forest loss dominates habitat loss
  hloss_full[inds[which(rowSums(cbind(nat2_full, nat1_full)) < 0 & nat1_full < nat2_full)]] <- 1
  
  #Other habitat loss dominates habitat loss
  hloss_full[inds[which(rowSums(cbind(nat2_full, nat1_full)) < 0 & nat2_full < nat1_full)]] <- 2
  
  #OBSERVED
  #Forest loss dominates habitat loss
  hloss_obs[inds[which(rowSums(cbind(nat2_full, nat1_full)) < 0 & nat1_obs < nat2_obs)]] <- 1
  
  #Other habitat loss dominates habitat loss
  hloss_obs[inds[which(rowSums(cbind(nat2_full, nat1_full)) < 0 & nat2_obs < nat1_obs)]] <- 2
  
  agrexp_hloss[[i]] <- list(stack(agrexp_obs, agrexp_full), stack(pastureexp_obs, pastureexp_full))
}

agrexp_hloss2 <- list(list(), list())
for(i in 1:6){
  agrexp_hloss2[[1]][[i]] <- agrexp_hloss[[i]][[1]]
  agrexp_hloss2[[2]][[i]] <- agrexp_hloss[[i]][[2]]
}

diff_list <- list()
figs4 <- list()
j <- 1
for(j in 1:2){
  cs_list <- agrexp_hloss2[[j]]
  diff <- matrix(rep(NA, 30), ncol = 5, nrow = 6)
  for(i in 1:length(cs_list)){
    ct <- crosstabm(cs_list[[i]][[1]], cs_list[[i]][[2]], percent = TRUE)
    diff[i,1] <- overallExchangeD(ct)
    diff[i,2] <- overallShiftD(ct)
    diff[i,3] <- overallQtyD(ct)
    diff[i,4] <- overallAllocD(ct)
    diff[i,5] <- overallDiff(ct)
  }
  diff_metrics <- as.data.frame(diff)
  colnames(diff_metrics) <- c("Exchange", "Shift", "Quantity", "Allocation", "Overall")
  metrics <- c("Allocation", "Quantity")
  diff_metrics$timestep <- ts
  diff_metrics$xpos <- as.factor(rownames(diff_metrics))
  diff_metrics %>%
    gather(value = "value", key = "difference", metrics) %>%
    mutate(difference = factor(difference, levels = metrics)) %>%
    ggplot(aes(x = xpos, y = value)) +
    geom_bar(. %>% filter(difference %in% metrics), stat = "identity", mapping = aes(fill = difference)) +
    scale_fill_manual(values = plasma(3, begin = 0.2, end = 0.8, direction = 1)) +
    scale_x_discrete(waiver(), labels= as.character(ts)) +
    scale_y_continuous(position = "left", limits = c(0,22)) + 
    xlab("Year") +
    ylab("Difference [%]") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = if(j==1) {element_blank()} else if(j==2) element_text(),
      axis.title.x = if(j==1) {element_blank()} else if(j==2) element_text(),
      axis.title.y = if(j==1) {element_text(colour = "transparent")} else if(j==2) {element_text(hjust = 1.5)},
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent",colour = NA),
      legend.position = if(j==1) {"none"} else if(j==2) c(0.7, 0.9),
      legend.direction="vertical",
      plot.margin = margin(t=0, r=10, b=10, l=10)
    ) -> 
    figs4[[j]]
}


boundary <- readOGR("/Users/simon/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/Amazon boundaries/Lim_Biogeografico.shp")
boundary <- gSimplify(boundary, tol = 0.05)

#Make levelplots of change rasters
s <- stack(agrexp_hloss[[6]])

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
                   axis.line=list(col='transparent')),
                 scales=list(draw=FALSE),            
                 col.regions=plasma(3, begin = 0.2, end = 0.8, direction = -1),
                 names.attr=rep('', nlayers(r))) +           
    layer(sp.polygons(boundary, lwd=2))
  
  if(i%in%c(1,3)){
    l <- l + layer(panel.text(-56, 9, labels = "observed", cex = 0.8, adj = 0, font = 2))
  }
  if(i%in%c(2,4)){
    l <- l + layer(panel.text(-56, 9, labels = "predicted", cex = 0.8, adj = 0, font = 2))
  }
  maps[[i]] <- l
}

#Make legends to be plotted seperately
legends_fig3 <- list()
for(i in 1:2){
  
  if(i == 1){
    lab <- c("no increase", "causing forest loss", "causing other habitat loss")
  }
  if(i == 2){
    lab <- c("no increase", "causing forest loss", "causing other habitat loss")
  }
  
  leg_data <- tibble(label = lab)
  leg_data$label <- factor(leg_data$label, lab)
  
  my_hist <- ggplot(leg_data, aes(label, fill = label)) + 
    geom_bar() + 
    scale_fill_manual(values= plasma(3, begin = 0.2, end = 0.8, direction = -1)) + 
    theme(title = element_blank(), 
          legend.position="bottom", 
          legend.margin=margin(t=-0, r=-0, b=10, l=-0, unit="pt"),
          legend.text = element_text(size = 8, colour = "black"))
  legends_fig3[[i]] <- cowplot::get_legend(my_hist)
}

png(file.path(figure_path, "figure4.png"), height = 140, width = 200, units = "mm", res = 1000)
plot_grid(
  plot_grid(figs4[[1]], figs4[[2]], ncol = 1, labels= c("a)", "b)"), hjust = c(-1.5,-1.5), vjust = c(1.5, 0.5)),
  plot_grid(
    plot_grid(maps[[1]], maps[[2]], ncol = 2), 
    legends_fig3[[1]], 
    plot_grid(maps[[3]], maps[[4]], ncol = 2), 
    legends_fig3[[2]], nrow = 4, 
    
    rel_heights = c(1/1.25, 1/12, 1/1.25, 1/12), 
    labels = c("c)", "", "d)", ""), vjust = c(1.5,1.5, 0.5, 1.5)), 
  rel_widths = c(1,3)
)
dev.off()

r <- mask
diff2 <- diff1 <- r
diff1[inds] <- rowSums(cbind(crop1_obs, crop2_obs, pasture_obs))
diff2[inds] <- rowSums(cbind(nat1_obs, nat2_obs))
diffstack <- stack(stack(diff2, diff1))
mind <- min(getValues(diffstack)[inds,])
maxd <- max(getValues(diffstack)[inds,])
maps2 <- list()
plot_seq <- seq(-1, 1, by = 0.125)

for(i in 1:nlayers(diffstack)){
  l <- levelplot(diffstack[[i]],
                 at=plot_seq,
                 margin=FALSE,                       
                 colorkey=FALSE,
                 xlab="", ylab="",
                 par.settings=list(
                   strip.border=list(col='transparent'),
                   panel.background = list(col = "transparent"),
                   strip.background=list(col='transparent'),
                   axis.line=list(col='transparent')
                 ),
                 scales=list(draw=FALSE),
                 col.regions= plasma(30, begin = 0.2, end = 0.8)) +
    layer(sp.polygons(boundary, lwd=2))
  if(i == 2){
    l <- l + layer(panel.text(-58, 9, labels = "Agr. expansion", cex = 0.7, adj = 0, font = 2))
  }
  if(i == 1){
    l <- l + layer(panel.text(-58, 9, labels = "Habitat loss", cex = 0.7, adj = 0, font = 2))
  }
  maps2[[i]] <- l
}

leg_data <- tibble(y = 1:length(plot_seq), change = plot_seq)
my_hist <- ggplot(leg_data, aes(x = change * 100, y = y, fill = change * 100)) + 
  geom_point() +
  guides(fill = guide_colorbar(title.position = "top",
                               title.hjust = 0.5,
                               barwidth = 8,
                               barheight = 1)) +
  scale_fill_viridis(option = "plasma", begin = 0.2, end = 0.8) +
  labs(fill = "Change [%]") +
  theme(legend.position="bottom",
        legend.box="horizontal",
        legend.text = element_text(size = 8, colour = "black"),
        legend.margin=margin(t=0, r=3, b=0, l=-3, unit="cm"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))

legend_fig2 <- cowplot::get_legend(my_hist)

sa <- readOGR("/Users/simon/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/countries_shp/countries.shp")
sa <- sa[sa$CONTINENT == "South America",]
sa <- gSimplify(sa, tol = 0.2)
sa <- st_as_sf(sa)
boundary_sf <- st_as_sf(boundary)
cent <- st_centroid(boundary_sf)
annotation <- data.frame(st_coordinates(st_cast(cent)), "label" = "Amazon catchment")

fig2a <- ggplot(data = sa) + 
  geom_sf(data = boundary_sf, fill = "darkgrey", colour = "darkgrey") +
  geom_sf(color = "black", fill = "transparent") +
  geom_text(data = annotation, aes(x = X, y = Y, label = label, fontface = 2), size = 3.5, hjust = 0.2, vjust = 1) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  coord_sf()

demands <- readRDS("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter2/data/data_ama/demands.rds")
dmd <- demands[c(1,27),-1]
length(which(!is.na(mask[])))
nrow(dat)
# r[r%in%crop] <- 1
# r[r%in%crop_mosaic] <- 2
# r[r%in%forest] <- 3
# r[r%in%grass] <- 4
# r[r%in%shrub] <- 5
# r[r%in%wetland] <- 6
# r[r%in%urban] <- 7
# r[r%in%other] <- 8
# r[r%in%water] <- 9

dmd_plotdata <- cbind(rowSums(dmd[,c(1:2)]),rowSums(dmd[,c(4:5)]), dmd[,3], rowSums(dmd[,c(6,8)]))
dmd_plotdata <- as.data.frame(rbind(dmd_plotdata, diff(dmd_plotdata))) * 100
classes <- c("Cropland",  "Pasture (grass/shrub)", "Forest", "Other habitat")
colnames(dmd_plotdata) <- classes

dmd_annotate <- round(dmd_plotdata[-3,], 3)
dmd_annotate[1,] <- paste0(dmd_annotate[1,], "%")
dmd_annotate[2,] <- paste0(dmd_annotate[2,], "%")

dmd_plotdata <- cbind(dmd_plotdata, "type" = c("1992", "2018", "difference"))
fig2b <- dmd_plotdata %>% 
  gather(key = "class", value = "change", -ncol(.)) %>% 
  mutate(class = factor(class, levels = classes)) %>%
  arrange(type) %>%
  mutate(v1992 = paste(round(change[which(type==1992)],2), "%")) %>%
  mutate(v2018 = paste(round(change[which(type==2018)],2), "%")) %>%
  filter(type == "difference") %>%
  mutate(y1992 = seq(1.25, (1.25 + nrow(.)-1), by = 1)) %>%
  mutate(y2018 = seq(0.75, (0.75 + nrow(.)-1), by = 1)) %>%
  
  ggplot(aes(y = class, x = change), fill = "black") + 
  geom_bar(stat = "identity") +
  xlim(-5.7, 5.7) +
  xlab("Change [%]") +
  ylab(NULL) +
  geom_text(aes(x = 0, y = y1992, label=v1992, hjust = ifelse(change > 0, 1,0)), size = 2.5, position = position_nudge(y = -0.05, x = 0.1 * c(-1, -1, 1, -1, -1, -1, -1))) +
  geom_text(aes(x = 0, y = y2018, label=v2018, hjust = ifelse(change > 0, 1,0)), size = 2.5, position = position_nudge(y = 0.05, x = 0.1 * c(-1, -1, 1, -1, -1, -1, -1))) +
  geom_text(x = 6, y = 7.2, label = "top: 1992", size = 2.5, hjust = 1) +
  geom_text(x = 6, y = 6.85, label = "bottom: 2018", size = 2.5, hjust = 1) +
  
  theme_bw() +   
  
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1,1,0.5,0), "cm"),
    axis.ticks = element_blank()
  )



fig2bleg <-plot_grid(fig2b, legend_fig2, ncol = 1, rel_heights = c(1, 0.1))
fig2ab <- plot_grid(fig2a, fig2bleg, ncol = 2, labels = c("a)", "b)"))
fig2c <- plot_grid(maps2[[2]], maps2[[1]],  ncol = 2)

png(file.path(figure_path, "figure2.png"), height = 140, width = 140, units = "mm", res = 1000)
plot_grid(fig2ab, fig2c, nrow = 2, labels = c("", "c)"), rel_heights = c(2, 1.5))
dev.off()



#The observed chanages are at similar magnitude to what we see in CGE models for Aus, Vietnam


#-------------------#
#### DESCRIPTION ####
#-------------------#

# This code provides the full validation of our model and figures (except for Fig 1):
# - Validation of fractional allocations and direction of change (Fig. 3 a, b)
# - Spatial-blocks cross-validation of suitability model at 1km and 10km resolution (Fig. 3 c)
# - Validation of agricultural expansion on natural habitat (Fig. 4)
# - Map of study area and historic changes in selected land use classes through time (Fig. 2)

rm(list = ls())

#-----------------------------#
#### 1. Packages and paths ####
#-----------------------------#


#Packages
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
require(flutes)
require(extraDistr)
require(wrswoR)

# Paths to results and required data

# Source helper functions
source("./R/functions.R")

# File paths
figure_path <- "/Users/simon/OneDrive - The University of Melbourne/PhD/writing/papers/2020_lu_model/GEB/figures"
data_path <- file.path(getwd(), "data", "data_ama")
lu_all <- readRDS(file.path(data_path, "lu.rds")) # Observed land use data (for validation)
mask <- readRDS(file.path(data_path, "mask_ama.rds")) # Amazon mask

boundary <- readOGR("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter2/data/bnd_ama/Lim_Biogeografico.shp") # Study aarea boundary feature

preds_full <- readRDS(file.path("outputs", "preds_full.rds")) # full model results
preds_semi <- readRDS(file.path("outputs", "preds_semi.rds")) # semi-naive model results
preds_naive <- readRDS(file.path("outputs", "preds_naive.rds")) # naaive model results

# sizing

# GEB
single_width <- 3.11 #(79mm, single column)
medium_width <- 4.33 #110mm, 1.5 x page)
full_width <- 6.61 #168 mm, two column width

#-----------------#
#### 2. Fig. 3 ####
#-----------------#

# Fig. 3 a
# Reorganise observed land use data

K <- 9 # number of land use classes
ts <- c(5, 10, 15, 20, 25, 27) # prediction/valdiation time steps

lu_obs <- list()
start <- seq(1,ncol(lu_all), by = K)

for(i in 1:length(start)){
  lu_obs[[i]] <- lu_all[,start[i]:(start[i]+(K-1))]
  colnames(lu_obs[[i]]) <- paste0("lu", 1:K, ".obs")
  print(i)
}

q1 <- lu_obs[[1]] # initial observed time step
q <- lu_obs[-1] # other observed time steps

# Calculate RMSE/MAE

# Get the maximum difference on any cell between first and consequtive time steps (for grouping the figure into 4 panels
q.diff <- lapply(q, FUN = function(x){rowMaxs(as.matrix(abs(x-q1)))})

# Null model: first time step becomes "prediction" (null model assumes no change through time) 
p_sn <- q1
null_rmse <- lapply(X = q, FUN = function(x) sqrt(rowMeans((p_sn - x)^2)))
null_mae <- lapply(X = q, FUN = function(x) rowSums(abs(p_sn - x)))

# Full, semi, naive model 

full_rmse <- list()
semi_rmse <- list()
naive_rmse <- list()

full_mae <- list()
semi_mae <- list()
naive_mae <- list()

for(i in 1:length(q)){
  full_rmse[[i]] <- sqrt(rowMeans((preds_full[[i]] - q[[i]])^2))
  semi_rmse[[i]] <- sqrt(rowMeans((preds_semi[[i]] - q[[i]])^2))
  naive_rmse[[i]] <- sqrt(rowMeans((preds_naive[[i]] - q[[i]])^2))
  full_mae[[i]] <- sqrt(rowSums((preds_full[[i]] - q[[i]])^2))
  semi_mae[[i]] <- sqrt(rowSums((preds_semi[[i]] - q[[i]])^2))
  naive_mae[[i]] <- sqrt(rowSums((preds_naive[[i]] - q[[i]])^2))
  print(i)
}

# Raw results data frame for plotting
data_f1 <- cbind(tibble(
  "null_rmse" = null_rmse %>% unlist(), 
  "full_rmse" = full_rmse %>% unlist(), 
  "semi_rmse" = semi_rmse %>% unlist(), 
  "naive_rmse" = naive_rmse %>% unlist(), 
  "full" = full_rmse - null_rmse,
  "semi" = semi_rmse - null_rmse,
  "naive" = naive_rmse - null_rmse,
  "Year" = rep(ts, each = nrow(q1))),
  "max.diff" = q.diff %>% unlist()) %>%
  mutate("category" = cut(max.diff, breaks=c(-Inf, 0.005, 0.3, 0.5, Inf), 
                          labels=c("less than 0.5% change", "0.5% - 30% change", "30% - 50% change", "more than 50% change")))

# Summarize raw results
data_f1 <- data_f1 %>% gather(key = "model", value = "diffs", "full", "semi", "naive")
data_aggr <- summarySE(data_f1, measurevar = "diffs", groupvars = c("Year", "category", "model"), conf.interval=0.95)
data_aggr$model <- factor(data_aggr$model, c("full", "semi", "naive"))

# Calculate fraction of cells in smallest category for results section
data_aggr %>% filter(category == "less than 0.5% change", model == "full") %>% dplyr::select(N) %>% sum -> cat0
data_aggr %>% filter(model == "full") %>% dplyr::select(N) %>% sum -> catall
cat0/catall * 100

data_aggr[,c(5:8)] <- data_aggr[,c(5:8)] * 100

# Make plot
fig3a <- ggplot(data_aggr, aes(x=Year, y=diffs, colour=model)) + scale_colour_viridis_d(option = "plasma", begin = 0.2, end = 0.6)
fig3a <- fig3a + geom_point(shape = 16, size = 1.5) + 
  geom_line(size=0.5) + 
  facet_wrap(~category) + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  theme_bw() + 
  scale_x_continuous(breaks = seq(5, 27, by = 5)) + 
  scale_y_continuous(limits=c(-0.07 * 100, 0.03 * 100)) +
  ylab("RMSE diff x 100") + 
  xlab("timestep") + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill="transparent", colour = "transparent"),
        strip.text = element_text(size=8, colour="black"), 
        legend.position = c(.5, .65), 
        legend.direction = "horizontal",
        legend.text = element_text(size = 8, colour = "black"),
        legend.title =  element_blank())


# Fig 3. b

# When observed change is always 0, the true negative is 1 and true positive 0
preds_null <- list()
for(i in 1:length(preds_full)){
  preds_null[[i]] <- p_sn
}

# Calculate difference metrics (Pontius et al.)
aucs_null <- diff_metrics(lu_obs, preds_null, mask, percent = TRUE, reference = 1)
aucs_naive <- diff_metrics(lu_obs, preds_naive, mask, percent = TRUE, reference = 1)
aucs_semi <- diff_metrics(lu_obs, preds_semi, mask, percent = TRUE, reference = 1)
aucs_full <- diff_metrics(lu_obs, preds_full, mask, percent = TRUE, reference = 1)
aucs <- list(aucs_null, aucs_naive, aucs_semi, aucs_full)

# Convert list with diff metrics to data frame and add descriptive columns
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

# Make plot
fig3b <- ggplot(t1, aes(x = timestep, y = diff, colour = model)) +
  geom_point(position=position_dodge(0.5), shape = 16, size = 1.5) + 
  scale_colour_viridis_d(option = "plasma", begin = 0.2, end = 0.8) + 
  theme_bw() + 
  ylab("disagreement [%]") + 
  xlab("timestep") +
  theme(
    legend.position = c(0.4, -0.5),
    legend.direction = "horizontal",
    legend.title =  element_blank(),
    legend.text = element_text(size = 8, colour = "black"),
    legend.margin=margin(c(0,0,0,0)),
    plot.margin=unit(c(5.5, 5.5, 18, 5.5), 'points')
  )

# Fig 3. c

# Load data (10 km and also 1 km data)
data_paths <- c(file.path(getwd(), "data", "data_ama"), file.path(getwd(), "data", "data_ama_1k"))
validation_results <- list()

# loop through 1k and 10k landscapes
for(z in 1:2){
  
  # load data and variables
  data_path <- data_paths[[z]]
  dat <- readRDS(file.path(data_path, "cov.rds"))
  lu_all <- readRDS(file.path(data_path, "lu.rds"))
  mask <- readRDS(file.path(data_path, "mask_ama.rds")) #country mask
  inds <- which(!is.na(mask[]))
  
  ts <- 1991 + c(1, 5, 10, 15, 20, 25, 27)
  K <- ncol(lu_all)/length(ts)
  n <- nrow(lu_all)
  lu <- lu_all[,1:K]
  rm(lu_all)
  
  # calcualte neighbourhood covariates
  weights <- list(matrix(1/9, 3, 3, byrow= TRUE)) #size of window
  weights <- rep(weights, length.out = K)
  ln <- neighbourhood(lu, c(1:K), weights, mask, enr = TRUE)
  ln <- scale(ln)
  
  cent <- attr(ln,c("scaled:center"))
  scal <- attr(ln,c("scaled:scale"))
  
  # Correlation analysis and susbetting
  pa <- dat[,which(grepl("PA", colnames(dat)))]
  pa <- pa[,1]
  dat <- dat[,-which(grepl("PA", colnames(dat)))]
  dat <- scale(dat)
  data <- cbind(dat, ln)
  
  corre <- correlations(data, sub = nrow(data)/50)
  preds <- colnames(corre)
  data <- data[,colnames(data)%in%preds]
  preds <- colnames(data)
  
  # create spatial blocks raster (9 blocks)
  tras <- raster(crs = crs(mask), ext = extent(mask), nrow = 3, ncol = 3)
  tras[] <- 1:length(tras)
  tras <- projectRaster(tras, mask, method = "ngb")
  tras[is.na(mask)] <- NA
  
  inds <- which(!is.na(mask[]))
  tras <- tras[inds]
  
  folds <- sort(na.omit(unique(tras)))
  
  # Suitability model cross validation
  tsub <- sample(1:nrow(data), round(nrow(data)/50))
  
  data_subs <- data[tsub,]
  lu_subs <- lu[tsub,]
  tras_subs <- tras[tsub]
  
  # Build models 
  
  # set up model formulas for different models (env, neigh, env+neigh)
  partype <- c("env", "neigh", "both")
  
  # Set up model formulas for different models 
  forms <- list(paste(preds[-which(grepl(preds, pattern = "ts"))], collapse = "+"),
                paste(preds[which(grepl(preds, pattern = "ts"))], collapse = "+"),
                paste(preds, collapse = "+"))
  
  mod <- list()
  rmse <- list()
  
  # Loop through models (env, neigh, both)
  for (j in 1:length(partype)){
    rmse.mod <- data.frame("rmse" = 1:length(inds), "fold" = NA)
    coefs.mod <- list()
    
    # Loop through 9 spatial blocks (cv-folds)
    for(i in folds){
      test <- i
      test_inds <- which(tras%in%test)
      
      # train model on subset of whole study area excluding the current fold
      train <- folds[!folds%in%test]
      train_inds <- which(tras_subs%in%train)
      suitmod <- suitmodel(form = forms[[j]], lu = lu_subs[train_inds,], data = data_subs[train_inds,], resolution = 1000, model = FALSE, maxit = 1000, decay = 0)
      
      # predict model to current fold
      pred <- predict(suitmod, newdata = data[test_inds,], type = "probs")
      
      # caclulate rmse on each fold (we want to get the RMSE spread across fold means, not across individual cell values)
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

# Prepare validation results
validation_results <- readRDS("validation_results.rds")
rmse <- lapply(validation_results, FUN = function(x) {x[[1]]})
resos <- c("10", "1")
df_final <- list()

for(i in 1:2){
  dfout <- data.frame("env" = rmse[[i]][[1]][,1]-rmse[[i]][[3]][,1], 
                      "neigh" = rmse[[i]][[2]][,1]-rmse[[i]][[3]][,1],
                      "folds" = rmse[[i]][[1]][,2],
                      "resolution" = rep(resos[i], length(rmse[[i]][[1]][,1]))
  )
  df_final[[i]] <- dfout
}

df_final <- tibble(do.call("rbind", df_final))

# Make data frame for plotting
df_plot <- 
  df_final %>% 
  gather(key = model, value = rmse, "env", "neigh") %>%
  group_by(resolution, model, folds) %>% #calculate fold-wise rmse means for each res and model
  summarize_at("rmse", list(fold.means = mean)) %>%
  group_by(resolution, model) %>% #calculate mean of fold means and min/max range
  summarize_at("fold.means", list(mean = mean, min = min, max = max, sd = sd))

# Make plot
fig3c <- ggplot(df_plot, aes(x = resolution, y = mean*100, col = model)) + 
  geom_point(position=position_dodge(0.5), shape = 16, size = 1.5) + 
  scale_colour_viridis_d(option = "viridis", begin = 0.2, end = 0.8) + 
  geom_errorbar(aes(ymin=(mean-sd)*100, ymax=(mean+sd)*100), width=.1, position=position_dodge(0.5)) + 
  theme_bw() +
  ylab("RMSE diff x 100") + 
  xlab("Resolution [km^2]") +
  theme(
    legend.position = c(0.4, -0.5),
    legend.direction = "horizontal",
    legend.title =  element_blank(),
    legend.text = element_text(size = 8, colour = "black"),
    legend.margin=margin(c(0,0,0,0)),
    plot.margin=unit(c(5.5, 5.5, 18, 5.5), 'points')
  )

# Combine panels into single figure
fig3bc <- plot_grid(fig3b, fig3c, nrow = 1, labels = c("(b)", "(c)"), vjust = c(0.2, 0.2), hjust = c(-1.1, -1.1), label_size = 12, rel_widths = c(2.5, 2))

pdf(file.path(figure_path, "figure_3.pdf"), height = 4.5, width = medium_width)
plot_grid(fig3a, fig3bc, nrow = 2, labels = c("(a)", ""), vjust = c(1.5, -1), rel_heights = c(3,2), label_size = 12)
dev.off()

#--------------#
#### Fig. 4 ####
#--------------#

# Determine agricultural land (pasture, crop) increase on natural habitat types (forest, grass/shrub etc) (both versions)
inds <- which(!is.na(mask[]))
agrexp_hloss <- list()

# Loop through predicted time steps
for (i in 1:length(preds_full)){
  
  r <- mask
  r[inds] <- 0
  
  # Estimate observed and predicted fractional land use change between time steps
  diff_obs <- lu_obs[[i+1]] - lu_obs[[1]] #observed
  diff_full <- preds_full[[i]] - lu_obs[[1]] #predicted
  
  # Identify observed fractional chaanges in agricultural classes
  crop1_obs <- diff_obs[,c(1)] #cropland 1
  crop2_obs <- diff_obs[,c(2)] #cropland 2
  pasture_obs <- diff_obs[,c(4,5)] #pasture
  
  # Identify observed fractional changes in natural habtiat classes
  nat1_obs <- diff_obs[,c(3)] # forest
  nat2_obs <- rowSums(diff_obs[,c(6,8)]) # naturaal habitat
  
  # Identify predicted fractional changes in agricultural classes
  crop1_full <- diff_full[,c(1)] #cropland 1
  crop2_full <- diff_full[,c(2)] #cropland 2
  pasture_full <- diff_full[,c(4,5)] #pasture
  
  # Identify predicted fractional changes in natural habitat classes
  nat1_full <- diff_full[,c(3)] # forest
  nat2_full <- rowSums(diff_full[,c(6,8)]) # natural habitat
  
  pastureexp_obs <- pastureexp_full <- agrexp_full <- agrexp_obs <- r
  
  # Reclassify to three classes: no change, agricultural expansion into forest, agricultural expansion into other natural hab types. We do this twice, for observed and for predicted time series.
  
  # Observed
  # Pasture cuts into forest
  pastureexp_obs[inds[which(pasture_obs > 0 & rowSums(cbind(nat1_obs,nat2_obs)) < 0 & nat1_obs < nat2_obs)]] <- 1
  
  # Pasture cuts into other hab types
  pastureexp_obs[inds[which(pasture_obs > 0 & rowSums(cbind(nat1_obs,nat2_obs)) < 0 & nat2_obs < nat1_obs)]] <- 2
  
  # Predicted
  # Pasture cuts into forest
  pastureexp_full[inds[which(pasture_full > 0 & rowSums(cbind(nat1_full,nat2_full)) < 0 & nat1_full < nat2_full)]] <- 1
  
  # Pasture cuts into other hab types
  pastureexp_full[inds[which(pasture_full > 0 & rowSums(cbind(nat1_full,nat2_full)) < 0 & nat2_full < nat1_full)]] <- 2
  
  # Observed
  # Cropland cuts into forest
  agrexp_obs[inds[which(rowSums(cbind(crop1_obs, crop2_obs)) > 0 & rowSums(cbind(nat1_obs,nat2_obs)) < 0 & nat1_obs < nat2_obs)]] <- 1
  
  # Cropland cuts into other hab types
  agrexp_obs[inds[which(rowSums(cbind(crop1_obs, crop2_obs)) > 0 & rowSums(cbind(nat1_obs,nat2_obs)) < 0 & nat2_obs < nat1_obs)]] <- 2
  
  # Predicted
  # Cropland cuts into forest
  agrexp_full[inds[which(rowSums(cbind(crop1_full, crop2_full)) > 0 & rowSums(cbind(nat1_full,nat2_full)) < 0 & nat1_full < nat2_full)]] <- 1
  
  # Cropland cuts into other hab types
  agrexp_full[inds[which(rowSums(cbind(crop1_full, crop2_full)) > 0 & rowSums(cbind(nat1_full,nat2_full)) < 0 & nat2_full < nat1_full)]] <- 2
  
  agrexp_hloss[[i]] <- list(stack(agrexp_obs, agrexp_full), stack(pastureexp_obs, pastureexp_full))
}

names(agrexp_hloss) <- c("agr_obs", "agr_full", "pas_obs", "pas_full")

#colorpal <- colorspace::diverging_hcl(4, h = c(250, 10), c = 70, l = c(26, 90), power = c(0.7, 1.7))
colorpal <- c("lightgrey", "darkgrey", plasma(2, begin = 0.4, end = 0.6, direction = -1))
colorpal <- c("lightgrey", "darkgrey", "orange3", "darkred")

# Find where observed and predicted maps of agricultural expansion agree and disagree (cor rej, false alarms, wrong hits, misses()
agrexp <- lapply(agrexp_hloss, FUN = function(x) {stack(x)})
agrexp_rasters <- list()
for(j in 1:6){
  overlays <- list()
  s <- agrexp[[j]]
  for(i in c(1,3)){
    overl <- mask
    overl[] <- NA
    overl[which(s[[i]][] == 0 & s[[i+1]][] == 0)] <- 1 #cor rejection
    overl[which(s[[i]][] == 1 & s[[i+1]][] == 1)] <- 2 #hit
    overl[which(s[[i]][] == 2 & s[[i+1]][] == 2)] <- 2 #hit
    overl[which(s[[i]][] == 0 & s[[i+1]][] == 1)] <- 3 #false alarm
    overl[which(s[[i]][] == 0 & s[[i+1]][] == 2)] <- 3 #false alarm
    overl[which(s[[i]][] == 1 & s[[i+1]][] == 2)] <- 4 #wrong hit (so rare that we lump them with misses)
    overl[which(s[[i]][] == 2 & s[[i+1]][] == 1)] <- 4 #wrong hit (so rare that we lump them with misses)
    overl[which(s[[i]][] == 1 & s[[i+1]][] == 0)] <- 4 #miss
    overl[which(s[[i]][] == 2 & s[[i+1]][] == 0)] <- 4 #miss
    overlays <- c(overlays, overl)
  }
  s_out <- stack(overlays)
  names(s_out) <- paste0("y_", ts[j], c("_agr", "_pas"))
  agrexp_rasters[[j]] <- s_out
}
agrexp_rasters <- stack(agrexp_rasters)

# Fig 4 a, b
exp_df <- agrexp_rasters[]
exp_df <- apply(na.omit(exp_df), 2, FUN = table)
exp_df <- exp_df/colSums(exp_df) * 100
exp_df <- exp_df %>%
  as.tibble %>%
  gather(key = "map", value = "frequency") %>% 
  separate(col = "map", sep = "_", into = c("y", "timestep", "type")) %>% 
  select(-"y") %>%
  add_column("error" = rep(as.factor(as.character(1:4)), 12))
exp_df$timestep <- factor(exp_df$timestep, levels = c(5, 10, 15, 20, 25, 27))
typs <- c("agr", "pas")
figs4 <- list()
for (j in 1:2){
  exp_df %>% 
    filter(type == typs[j]) %>% 
    ggplot(aes(x = timestep, y = frequency, fill = factor(error, levels = c("4", "3", "2", "1")))) + 
    geom_bar(stat = "identity") +
    scale_fill_manual(values = rev(colorpal)) +
    labs(fill = "") +
    ylab("[%] of landscape") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(),
      axis.text.y = element_text(margin = margin(t = 0, r = -1, b = 0, l = 0)),
      axis.title.x = element_text(),
      axis.title.y = element_text(hjust = 0.5, vjust = -1),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "none",
      legend.direction="vertical",
      plot.margin = margin(t=0, r=10, b=10, l=10)
    ) -> figs4[[j]]
}

# Fig 4 c, d
lattice.options(
  layout.heights=list(bottom.padding=list(x=-1), top.padding=list(x=0)),
  layout.widths=list(left.padding=list(x=-0.5), right.padding=list(x=9))
)

maps <- list()
for (i in 1:2){
  r <- agrexp_rasters[[c(11,12)]][[i]]
  r <- ratify(r)
  rat <- levels(r)[[1]]
  rat$ID <- c("1", "2", "3", "4")
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
                 col.regions = colorpal,
                 names.attr=rep('', nlayers(r))) +           
    layer(sp.polygons(boundary, lwd=1))
  maps[[i]] <- l
}


# Make legend
legends_fig4 <- list()

lab <- c("Correct: observed persistence predicted as persistence", 
         "Correct: predicted agricultural increase results in decrease \n in correct habitat class", 
         "Error: observed persistence predicted as agricultural increase", 
         "Error: observed agricultural increase predicted as persistence \n or results in decrease in incorrect  habitat class")

leg_data <- tibble(label = lab)
leg_data$label <- factor(leg_data$label, lab)

my_hist <- ggplot(leg_data, aes(label, fill = label)) + 
  guides(fill=guide_legend(ncol=1)) +
  geom_bar() + 
  scale_fill_manual(values= colorpal) + 
  theme(title = element_blank(), 
        legend.position="left", 
        legend.margin=margin(t=0, r=10, b=10, l=0, unit="pt"),
        legend.text = element_text(size = 8, colour = "black"))
legends_fig4 <- cowplot::get_legend(my_hist)

mp <- list()

# Constract sub panels
labels <- list(c("(a)", "(b)"), c("(c)", "(d)"))

for (i in 1:2){
  mp[[i]] <- ggdraw(maps[[i]]) +
    draw_plot(figs4[[i]], .62, 0, .38, .85) +
    draw_plot_label(
      labels[[i]],
      c(0, 0.65),
      c(1, 1),
      size = 12
    )
}

# Combine Fig 4. subpanels
pdf(file.path(figure_path, "figure_4.pdf"), height = 6, width = medium_width)

plot_grid(mp[[1]], mp[[2]], ncol = 1, legends_fig4, nrow = 3, rel_heights = c(1, 1, 0.5), labels = c("cropland", "pasture"), label_size = 12, vjust = c(18, 18), hjust = c(-3.3, -4))

dev.off()

#--------------#
#### Fig. 2 ####
#--------------#

# Fig 2. a

# South America
sa <- readOGR("/Users/simon/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/countries_shp/countries.shp")
sa <- sa[sa$CONTINENT == "South America",]
sa <- gSimplify(sa, tol = 0.2)
sa <- st_as_sf(sa)
boundary_sf <- st_as_sf(boundary)
cent <- st_centroid(boundary_sf)
annotation <- data.frame(st_coordinates(st_cast(cent)), "label" = "Amazon catchment")
annotation[,2] <- annotation[,2] - 0.5

# Add protected areas
pa <- raster("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter2/data/data_ama/temp/PA1992_ama.tif")
pa[which(pa[] == 1)] <- NA
pa <- gplot_data(pa)

# Make figure 2 a
fig2a <- ggplot(data = sa) + 
  geom_sf(data = boundary_sf, fill = "darkgrey", colour = "darkgrey") +
  geom_sf(color = "black", fill = "transparent") +
  geom_tile(data = dplyr::filter(pa, !is.na(value)), 
            aes(x = x, y = y), fill = "black") + 
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

# Extract by how much forest, pasture and cropping classes change (this is the aggregated change from the observed tmie series)
demands <- readRDS("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter2/data/data_ama/demands.rds")
dmd <- demands[c(1,27),-1]
dmd_plotdata <- cbind(rowSums(dmd[,c(1:2)]),rowSums(dmd[,c(4:5)]), dmd[,3])
dmd_plotdata <- as.data.frame(rbind(dmd_plotdata, diff(dmd_plotdata))) * 100
classes <- c("Cropping",  "Pasture", "Forest")
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
  theme_bw() +   
  
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1,0.5,1,0), "cm"),
    axis.ticks = element_blank()
  )


pdf(file.path(figure_path, "figure_2.pdf"), height = 3, width = single_width)
plot_grid(fig2a, fig2b, ncol = 2, labels = c("(a)", "(b)"), label_size = 12)
dev.off()
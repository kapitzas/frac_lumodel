####-------------------------------####
#### Some plotting/helper functions####
####-------------------------------####

# summarySE is from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  colnames(datac)[which(colnames(datac) == "mean")] <-  measurevar
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

# Calculate difference metrics according to Pontius 2011 (for Fig. 3. b) 
diff_metrics <- function(obs, preds, mask, reference = NULL,...){
  K <- ncol(lu_obs[[1]])
  preds <- c(obs[1], preds[1:6])
  cont_tables <- list()
  inds <- which(!is.na(mask[]))
  out <- list()
  for(j in 1:K){
    change_obs <- change_pred <- mask
    for (i in 2:(length(obs))){
      if(!is.null(reference)){
        ref <- reference
      }else{
        ref <- i-1
      }
      change_pred[inds[which(preds[[ref]][,j] > preds[[i]][,j])]] <- 1 #decrease
      change_pred[inds[which(preds[[ref]][,j] == preds[[i]][,j])]] <- 2 #same
      change_pred[inds[which(preds[[ref]][,j] < preds[[i]][,j])]] <- 3 #increase
      
      change_obs[inds[which(obs[[ref]][,j] > obs[[i]][,j])]] <- 1 #decrease
      change_obs[inds[which(obs[[ref]][,j] == obs[[i]][,j])]] <- 2 #same
      change_obs[inds[which(obs[[ref]][,j] < obs[[i]][,j])]] <- 3 #increase
      
      cont_tables[[i-1]] <- diffeR::crosstabm(change_pred, change_obs, ...)
    }
    out[[j]] <- cont_tables
  }
  out
}

# this function is from https://stackoverflow.com/questions/47116217/overlay-raster-layer-on-map-in-ggplot2-in-r
gplot_data <- function(x, maxpixels = 50000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- dplyr::as.tbl(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}


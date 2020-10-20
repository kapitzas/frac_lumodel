####--------------------------------####
####1) Calculate neighbourhood layer####
####--------------------------------####

neighbourhood <- function(lu, cols, weights, mask, ..., enr = F, suffix = NULL){
  if(length(which(!is.na(mask[]))) != nrow(lu)) stop("mask and land use not same length")
  out <- lu[,cols]
  for (i in 1:length(cols)){
    w <- weights[[cols[i]]]
    l <- mask
    inds <- which(!is.na(getValues(l)))
    l[inds] <- lu[,cols[i]]
    l <- raster::focal(l, w, na.rm = TRUE, pad = TRUE,...)
    if(enr == TRUE){ #Should ouput be converted to enrichment factors (according to Verburg et al 2004)
      l <- l/mean(lu[,cols[i]])
    }
    out[,i] <- getValues(l)[inds]
  }
  if(!is.null(suffix)){
    colnames(out) <- paste0("lu_", cols, "_", suffix)
  }
  out
}

####--------------------------------####
####2) Automated predictor selection####
####--------------------------------####

correlations <- function(covs, thresh = 0.7, sub = inds) {
  subs_cor <- sample(1:nrow(covs), size = sub) #sample subset to do corr analysis on
  cors <- cor(covs[subs_cor, ], method = "spearman") #get pred correlations based on that subset
  while (min(abs(cors[abs(cors) >= thresh])) != 1){ #reduce predictor set until only uncorrelated pairs remain
    values <- cors[which(abs(cors) > thresh)]
    corellated <- which(abs(cors) > thresh)
    values[values ==1] <- NA
    corellated[which(values== max(values, na.rm = T))]
    rows_highest_cor <- which(cors == max(values, na.rm = T), arr.ind = T)[,1]
    cors_cur <- abs(cors[rows_highest_cor,])
    '%ni%' <- Negate('%in%')
    m1 <- max(cors_cur[1,][cors_cur[1,]%ni%c(max(values, na.rm = T),1)])
    m2 <- max(cors_cur[2,][cors_cur[2,]%ni%c(max(values, na.rm = T),1)])
    out <- ifelse(m1 > m2, 1, 2)
    cors <- cors[-which(colnames(cors) == names(rows_highest_cor)[out]), -which(colnames(cors) == names(rows_highest_cor)[out])]
    nrow(cors)
  }
  return(cors)
}

####--------------------------####
####3) Build suitability model####
####--------------------------####

suitmodel <- function(form, lu, data, sub = NULL, resolution,...){
  
  #Sample random subset of data to build model on
  if(!is.null(sub)){
    subs <- sample(1:nrow(lu), sub)
  }else{
    subs <- 1:nrow(lu)
  }
  
  #Turn fractions into frequency counts for multinomial regression
  counts <- integerify(x = lu[subs,], resolution = resolution)
  data_sub <- as.data.frame(data[subs,])
  
  #Pass pre-determined formula
  f <- as.formula(paste("counts", "~", form))
  suit_model <- nnet::multinom(f, data = data_sub, ...)
  suit_model
}

####-----------------------------------------------####
####4) Extract demand change from observed land use####
####-----------------------------------------------####

demand <- function(inds = NULL, landuse, ts, path = NULL, k, type = "mean"){
  lu <- landuse
  if(!is.null(inds)){
    lu <- lu[inds, ] #get subset to extract demand from (i.e. when cells belong to specific regions/countries)
  }
  
  #Get the current land use supply in each class
  if(type == "mean"){
    class_supply <- as.numeric(colMeans(lu))
  }
  
  if(type == "sum"){
    class_supply <- colSums(lu)
  }
  
  years <- ts[1]:tail(ts,1)
  demand <- matrix(NA, nrow = length(years), ncol = k + 2)
  demand[,1] <- years
  obs_ind <- which(demand[,1]%in%ts)
  
  for(i in 1:(ncol(lu)/k)){
    demand[obs_ind[i], 2:(k+1)] <- class_supply[seq((i * k - (k-1)), i * k, by = 1)]
  }
  
  #interpolate intended time steps
  for(i in 1:k){
    demand[, i + 1] <- approx(demand[,1], demand[,i + 1], xout = years)$y
  }
  
  demand[, k + 2] <- rowSums(demand[,-1], na.rm = TRUE)
  if(!is.null(path)){
    saveRDS(demand, path)
  }
  return(demand)
}

####------------------####
####5) Allocate demand####
####------------------####
allocation <- function(lu, sm, params, dmd, ln, constraint, pa = NULL){
  
  #number of land use classes and number of cells
  k <- ncol(lu)
  n <- nrow(lu)
  
  resolution <- params$resolution
  max_dev <- params$max_dev
  growth <- params$growth
  no_change <-  1:K%in%params$no_change
  
  #Turn intiital land use into integer counts
  p_t0 <- integerify(lu, resolution = resolution)
  
  #Convert demand time series to integer counts
  supply_t0 <- colSums(p_t0) #Initial land use
  demand_traj <- integerify(x = dmd, resolution = sum(supply_t0))
  demand_t1 <- demand_traj[2,] #Second row of demand trajectory is demand to be allocated
  supply_t1_candidate <- demand_traj[1,] #First time step land use supply becomes "candidate" that gets recalculated in each iteration until it meets demand_t1
  
  dev_diff <- abs(diff(demand_traj))/demand_traj[1,] * 100 #Initital candidate % deviation of current supply from demand, gets recalculated in each iteration until it's below max_dev
  dev_diff[which(is.na(dev_diff))] <- 0
  
  #Counter
  count <- 0
  
  if(!is.null(pa)){
    pa_inds <- which(pa==0)
    supply_t0_pa <- colSums(p_t0[pa_inds,])
    nopa_inds <- which(pa==1)
    p_t0 <- p_t0[nopa_inds,]
    ln <- ln[nopa_inds,]
    sm <- sm[nopa_inds,]
  }
  
  #First land use map becomes "candidate" on which allocations take place iteratively
  p_t1_candidate <- p_t0
  
  #Determine which classes are 0 on each cell (so we can constrain growth to those that aren't 0)
  if(constraint){
    are_zero <- p_t0 == 0
    inds_list <- list()
    i <- 1
    for(i in 1:K){
      inds_zero_all <- which(are_zero[,i]) #all that are 0 should stay 0, except for...
      inds_zero <- inds_zero_all[which(ln[inds_zero_all,i]!=0)] # a subset of cells near where land use already exists (neigbourhood)
      size <- length(are_zero[,i]) * ((growth[i])/100) # we want to sample this many from that subset to remain.
      if(size < length(inds_zero)){
        keep <- inds_zero[sample_int_rank(length(inds_zero), size = size, prob = sm[inds_zero,i])]
        inds_list[[i]] <- inds_zero_all[-which(inds_zero_all%in%keep)]
      }
      if(size > length(inds_zero)){
        leftover <- size - length(inds_zero)
        inds_leftover <- which(!inds_zero_all%in%inds_zero)
        keep <- inds_leftover[-sample_int_rank(length(inds_leftover), size = min(c(leftover, length(inds_leftover))), prob = sm[inds_leftover,i])]
        inds_list[[i]] <- inds_zero_all[-which(inds_zero_all%in%c(inds_zero, keep))]
      }
    }
  }
  
  # Iterative allocation
  while (any(dev_diff > max_dev)) {
    
    #Counter increment
    count <- count + 1
    
    #Demand to be allocated
    demand_change <- demand_t1 - supply_t1_candidate
    
    if(!is.null(pa)){
      demand_change <- demand_change - supply_t0_pa
    }
    
    #Calculate change factor (a),by how much do we have to multiply the candidate land use proportions to satisfy the modelled suitability.
    ideal_change <- (sm * resolution) / p_t1_candidate
    
    both_0 <- is.na(ideal_change) #this determines which cells are 0 in the suitability map and the current iteration. We keep those as they are.
    ideal_change[both_0] <- 1
    
    cand_0 <- !is.finite(ideal_change) #This determines which cells are 0 in p_t1_candidate
    ideal_change[cand_0] <- sm[cand_0]
    
    #Calculate Relative suitability (r) from change factors. Sums to 1 in each column
    rel_suitability <- ideal_change %*% diag(1/colSums(ideal_change))
    
    #Allocate demand change between pixels (d), wieghted by r (d * r = m) and adjusted by stepsize (default is 1 but can be smaller for more fine-scale allocations: more stable, but slower)
    target_lu_change_pixel <-  rel_suitability %*% diag(demand_change)
    
    #Add changes to candidate map and make everyting positive.
    p_t1_proposal <- p_t1_candidate + target_lu_change_pixel
    p_t1_proposal <- pmax(p_t1_proposal, 0)
    
    #Where land use demand is 0, make all cells 0 in that land use
    if(any(demand_t1 == 0)){
      p_t1_proposal[, which(demand_t1 == 0)] <- 0
    }
    
    if(constraint){
      for(i in (1:K)){
        p_t1_proposal[inds_list[[i]], i] <- 0
      }
    }
    
    #Turn proposed land use map into integers
    p_t1_candidate <- integerify(x = p_t1_proposal,  resolution = resolution, no_decrease = no_change, z = p_t1_candidate)
    
    #ii Calcuate new candidate supply, i.e. the supply of the currently proposed candidate
    supply_t1_candidate <- colSums(p_t1_candidate)
    
    #Recalculate % deviation of candidate supply from demand
    diff <- abs(demand_t1 - supply_t1_candidate)
    if(!is.null(pa)){
      diff <- abs(demand_t1 - (supply_t1_candidate + supply_t0_pa))
    }
    dev_diff <- diff/demand_t1 * 100
    dev_diff[which(is.na(dev_diff))] <- 0
    
    cat("\r", paste0("Iteration: ", count, "    "), "Deviation from target per class [%]: ", paste(round(dev_diff, 3), sep = " "))
  }
  
  #When allocations are ready, return result
  if(!is.null(pa)){
    pred_out <- matrix(NA, ncol = K, nrow = n)
    pred_out[nopa_inds,] <- p_t1_candidate/resolution
    pred_out[pa_inds,] <- lu[pa_inds,]
  }else{
    pred_out <- p_t1_candidate/resolution
  }
  colnames(pred_out) <- colnames(lu)
  pred_out
}

####---------------------------------####
####6) Turn proportions into integers####
####---------------------------------####

#Integerify function
integerify <- function (x, z = NULL,  resolution = resolution, no_decrease = NULL) {
  
  n <- nrow(x)
  k <- ncol(x)
  
  if (any(no_decrease)) { #if there is at least one class that deosn't decrease (typically urban)
    min_allocation <- matrix(0, n, k)
    for (class in which(no_decrease)) { #loop through classes
      min_allocation[, class] <- z[, class]
    }
  } else {
    min_allocation <- NULL
  }
  
  # if a minimum number must be placed in a certain column, withold them and just allocate the others
  if (!is.null(min_allocation)) {
    to_allocate <- resolution - rowSums(min_allocation)
    
    # need to reduce the probability of allocating the remainder to these cells too
    # so make x sum to resolution, subtract the minima, and make sure the result is positive
    x <- x - min_allocation #might have to be other way round
    x <- resolution * x / rowSums(x) #* resolution because gets scaled to 0-1
    
    x <- pmax(x, 0)
    
  } else {
    to_allocate <- rep(resolution, n)
  }
  
  # random draws of the non-reserved ones
  x <- extraDistr::rmnom(n, to_allocate, prob = x)
  
  # add on the reserved ones again
  if(!is.null(min_allocation)) {
    x <- x + min_allocation
  }
  x
}


####---------------------------------####
####7) Some plotting/helper functions####
####---------------------------------####

makeRaster <- function(mask, lu, ts = NULL, class = NULL){
  if(!is.null(ts)){
    lu <- lu[[ts]]
  }
  if(!is.null(class)){
    mask[which(!is.na(mask[]))] <- lu[,class]
  }else{
    mask[which(!is.na(mask[]))] <- lu
  }
  
  mask
}

#Apply neighbourhood raster

elasticities <- function(change, elas){
  
  for(i in 1:length(elas)){
    fq <- (1-elas[i])/2
    sq <- elas[i] + fq
    qs <- quantile(change[,i], c(fq, sq))
    inds <- which(change[,i] >= qs[1] & change[,i] <= qs[2])
    change[inds,i] <- 0
  }
  change
}

elasticities2 <- function(ideal_change, elas){
  for(i in 1:length(elas)){
    size <- ceiling(nrow(ideal_change) * (1-elas[i]))
    inds <- sample(1:length(ideal_change[,i]), size = size)
    ideal_change[-inds,i] <- 0
  }
  ideal_change
}

#below is from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

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

#Caluclate metrics
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

#https://stackoverflow.com/questions/47116217/overlay-raster-layer-on-map-in-ggplot2-in-r
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


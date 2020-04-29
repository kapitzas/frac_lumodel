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
allocation <- function(lu, sm, params, dmd, ln, constraint){
  
  k <- ncol(lu)
  n <- nrow(lu)
  
  resolution <- params$resolution
  no_change <- (1:k)%in%params$no_change
  
  #Turn intiital land use into integer counts
  p_t0 <- integerify(lu, resolution = resolution)
  
  #Convert demand time series to integer counts
  supply_t0 <- colSums(p_t0) #Initial land use
  if(any(no_change)){
    sum_nochange <- sum(supply_t0) - supply_t0[which(no_change)]
    demand_traj <- integerify(x = dmd[,-which(no_change)], resolution = sum_nochange) #Turn into integer counts
    
    demand_traj <- cbind(demand_traj, supply_t0[which(no_change)])
  }else{
    demand_traj <- integerify(x = dmd, resolution = sum(supply_t0))
  }
  
  demand_t1 <- demand_traj[2,] #Second row of demand trajectory is demand to be allocated, first 
  
  #First time step land use supply becomes "candidate" that gets recalculated in each iteration until it meets demand_t1
  supply_t1_candidate <- supply_t0
  
  #number of land use classes and number of cells
  
  #Initital candidate % deviation of current supply from demand, gets recalculated in each iteration until it's below max_dev
  max_dev <- params$max_dev
  dev_diff <- abs(diff(dmd))/dmd[1,] * 100
  dev_diff[which(is.na(dev_diff))] <- 0
  stepsize <- params$stepsi
  ch_thresh <- params$ch_thresh
  
  #Counter
  count <- 0
  max_it <- params$max_it
  
  #First land use map becomes "candidate" on which allocations take place iteratively
  p_t1_candidate <- p_t0
  are_zero <- which(p_t0 == 0)
  
  # Iterative allocation
  while ((sum(dev_diff > rep(max_dev,k)) !=0) & count < max_it) {
    
    #Counter increment
    count <- count + 1
    
    #Demand to be allocated
    demand_change <- demand_t1 - supply_t1_candidate
    
    #Calculate change factor (a),by how much do we have to multiply the candidate land use proportions to satisfy the modelled suitability.
    ideal_change <- (sm * resolution) / p_t1_candidate
    
    both_0 <- is.na(ideal_change)
    ideal_change[both_0] <- 0
    
    cand_0 <- !is.finite(ideal_change)
    ideal_change[cand_0] <- sm[cand_0]
    
    ideal_change[which(ideal_change == 1)] <- 0
    
    #Calculate Relative suitability (r) from change factors. Sums to 1 in each column
    rel_suitability <- ideal_change %*% diag(1/colSums(ideal_change))
    
    #rel_suitability[,no_change] <- 0
    #Allocate demand change between pixels (d), wieghted by r (d * r = m) and adjusted by stepsize (default is 1 but can be smaller for more fine-scale allocations: more stable, but slower)
    target_lu_change_pixel <-  rel_suitability %*% diag(demand_change)  * stepsize
    
    #Where the current candidate deviation of supply from demand is smaller than maximum allowed deviation, make all changes 0.
    if(any(dev_diff <= max_dev)){
      mini_change  <- na.omit(dev_diff <= max_dev)
      target_lu_change_pixel[,mini_change] <- 0
    }
    
    #Add changes to candidate map and make everyting positive.
    p_t1_proposal <- p_t1_candidate + target_lu_change_pixel
    p_t1_proposal <- pmax(p_t1_proposal, 0)
    
    #Where land use demand is 0, make all cells 0 in that land use
    if(any(demand_t1 == 0)){
      p_t1_proposal[, which(demand_t1 == 0)] <- 0
    }

    if(!is.null(constraint)){
      if(constraint == "all"){
        p_t1_proposal[are_zero] <- 0
      }
    }
    
    #Turn proposed land use map into integers
    p_t1_candidate <- integerify(x = p_t1_proposal,  resolution = resolution, no_decrease = no_change, z = p_t1_candidate)
    
    #ii Calcuate new candidate supply, i.e. the supply of the currently proposed candidate
    supply_t1_candidate <- colSums(p_t1_candidate)
    
    #Recalculate % deviation of candidate supply from demand
    diff <- abs(demand_t1 - supply_t1_candidate)
    dev_diff <- diff/demand_t1 * 100
    dev_diff[which(is.na(dev_diff))] <- 0
    
    cat("\r", paste0("Iteration: ", count, "    "), "Deviation from target per class [%]: ", paste(round(dev_diff, 3), sep = " "))
    
    if (count == 20000) {
      stop("algorithm failed to converge")
    }
  }
  
  #When allocations are ready, return result
  pred_out <- p_t1_candidate/resolution
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

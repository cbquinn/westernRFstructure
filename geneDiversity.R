# function to calculate gene diversity for haplotype data according to Nei 1987 
# Nei M (1987) Molecular Evolutionary Genetics. Columbia University Press: New York

# x is a vector of haplotype IDs

# point estimate
geneDiversity <- function(x) {
  x <- x[!is.na(x)]
  unbiasN <- length(x)/(length(x)-1)
  k <- as.character(unique(x))
  squaredp <- rep(NA, length(k))
  names(squaredp) <- k
  for(i in k){
    squaredp[i] <- (length(x[x==i])/length(x))^2
  }
  unbiasN * (1-sum(squaredp)) 
}

# variance
geneDiversity_variance <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  k <- as.character(unique(x))
  
  term1 <- 2/(n*(n-1))
  term2 <- 2*(n-2)
  temp3 <- rep(NA, length(k))
  names(temp3) <- k
  for(i in k){
    temp3[i] <- (length(x[x==i])/length(x))^3
  }
  term3 <- sum(temp3)
  temp4 <- rep(NA, length(k))
  names(temp4) <- k
  for(i in k){
    temp4[i] <- (length(x[x==i])/length(x))^2
  }
  term4 <- sum(temp4)^2
  temp5 <- rep(NA, length(k))
  names(temp5) <- k
  for(i in k){
    temp5[i] <- (length(x[x==i])/length(x))^2
  }
  term5 <- sum(temp5)
  term1*((term2*(term3 - term4)) + term5 - term4)
}

# standard deviation
geneDiversity_sd <- function(x) {
  var <- geneDiversity_variance(x)
  var^(0.5)
}

# standard error
geneDiversity_se <- function(x) {
  var <- geneDiversity_variance(x)
  var^(0.5) / (length(x)^(0.5))
}

# function modified from sGD function in sGD package that calculates mitochondrial gene diversity by spatial neighborhoods 
# (original function is designed for diploid genotypes)

# mtHaplos -- dataframe with sampleID, and mitochondrial haplotype
# xy -- dataframe of coordinates
# dist.mat -- euclidian distance matrix
# NH_radius -- neighborhood radius
# min_N -- minimum number of individuals in neighborhood, otherwise value is returned as NA

mysGD_mtDNA <- function (mtHaplos, xy, dist.mat, NH_radius, min_N, max_N = NULL) {
  cat("Reading input files...\n")
  N <- nrow(mtHaplos)
  cat("Input summary:\n")
  cat(paste("\t individuals:", N, "\n"))
  cat(paste("\t minimum sample size:", min_N, "\n"))
  if (is.numeric(max_N)) {
    cat(paste("\t maximum sample size:", max_N, "\n"))
  }
  else {
    cat(paste("\t maximum sample size: NA", "\n"))
  }
  if (length(NH_radius) > 1) {
    min_R <- min(NH_radius, na.rm = T)
    max_R <- max(NH_radius, na.rm = T)
    cat(paste0("\t variable neighborhood radius (min = ", 
               min_R, ";max = ", max_R, ")", "\n"))
  }
  else {
    cat(paste("\t neighborhood radius:", NH_radius, "\n"))
  }
  cat("Determining neighborhood membership from dist.mat and NH_radius...\n")
  NHmat <- c()
  # make a matrix of neighborhood membership
  for (r in c(1:N)) {
    if (length(NH_radius) == 1) {
      radius <- NH_radius
    }
    else {
      radius <- NH_radius[r]
    }
    rvals <- ifelse(dist.mat[r, ] < radius, 1, 0) # for individual r, 1 means in the same neighborhood, 0 means not
    NHmat <- rbind(NHmat, rvals)
  }
  NH_N <- as.numeric(rowSums(NHmat, na.rm = T))   # neighborhood size for each individual
  if (max(NH_N) < min_N) {
    stop("There are no neighborhoods with at least min_N individuals at the radius selected. Aborting sGD.")
  }
  if (length(NH_radius) == 1) {
    valid_pops <- as.numeric(which(NH_N >= min_N))
  }
  else {
    valid_pops <- as.numeric(which(!is.na(NH_radius) & NH_N >= 
                                     min_N))
  }
  npops = length(valid_pops)
  NH_Index <- c(1:N)
  NH_summary <- data.frame(NH_Index, as.character(mtHaplos$sampleID), 
                           xy, NH_N, stringsAsFactors = F)
  names(NH_summary) <- c("NH_Index", "NH_ID", "X", "Y", "N")
  
  
  NH_summary$geneD <- NA
  # create empty list to fill u5
  for (pop in valid_pops) {
    if (!is.null(max_N) && NH_summary$N[pop] > max_N) {
      pop_ind <- sample(which(NHmat[pop, ] == 1), max_N)
    }
    else {
      pop_ind <- which(NHmat[pop, ] == 1)
    }
    haplos <- mtHaplos[pop_ind, 2]
    geneD <- geneDiversity(haplos)
    NH_summary$geneD[pop] <- geneDiversity(haplos)
  }
  NH_summary
}


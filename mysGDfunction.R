# see Shirk A, Cushman S (2011). sGD: software for estimating spatially explicit indices of genetic diversity. Molecular Ecology Resources 11(5): 922-934
# functions modified slightly to a) avoid errors thrown from a combination of loci with 2 and 3 alleles and b) to incorporate haplotype frequency data

# this takes a genind object in, and outputs a dataframe
# just like the original genind_obj, but deals better with the fact that some of my loci have 2 alleles and some have 3
genind2df_digcorrect <- function(genind_obj){
  df1 <-  genind2df(genind_obj, oneColPerAll = TRUE)
  #figure out which columns only have 2 digits (length of 2)
  genos <- df1 %>%
    dplyr::select(-pop)
  temp <- genos %>% 
    dplyr::mutate_all(as.numeric)
  
  twoers <- temp[temp < 100 & !(is.na(temp))]
  twoers_conv <- paste0("0",as.character(twoers))
  
  #replace
  genos[temp < 100 & !(is.na(temp))] <- twoers_conv
  
  loc1 <- seq(1, ncol(genos), 2)
  loc2 <- loc1+1
  
  newgenos <- as.data.frame(matrix(nrow=nrow(genos), ncol=31))
  
  for(i in 1:length(loc1)) {
    newgenos[,i] <- paste0(genos[,loc1[i]], genos[,loc2[i]])
    tempname <- names(genos)[loc1[i]]
    names(newgenos[,i]) <- gsub(pattern = "\\.1", "", tempname)
  }
  
  # replace NANA with NA
  newgenos[newgenos=="NANA"] <- NA
  
  #put sampleID row names and pop column back on
  df_dat <- cbind(pop=df1$pop, newgenos)
  rownames(df_dat) <- row.names(df1)
  return(df_dat)
}

# slight alterations to original function to deal with errors thrown from the fact that some of my loci have 2 alleles and some have 3 and some errors finding NeEstimator executable
mysGD <- function (genind_obj, xy, dist.mat, NH_radius, min_N, max_N = NULL, 
                 metrics = NULL, NHmat_ans = FALSE, genout_ans = FALSE, file_name = NULL, 
                 NeEstimator_dir = NULL) {
  if (packageVersion("adegenet") < "2.0.0") {
    stop("Please install the latest version of the adegenet package (>= 2.0.1)")
  }
  if (packageVersion("hierfstat") < "0.04.22") {
    stop("Please install the latest version of the hierfstat package (>= 0.04.22)")
  }
  cat("Reading input files...\n")
  df_dat <- genind2df_digcorrect(genind_obj)
  allele.digits <- nchar(genind_obj@all.names[[1]][1])
  df_dat[is.na(df_dat)] <- paste(rep("0", allele.digits * 2), 
                                 collapse = "")
  numloci <- length(names(genind_obj@all.names))
  N <- dim(genind_obj@tab)[1]
  OS <- as.character(Sys.info()["sysname"])
  genout_file <- paste(file_name, "_genepop.gen", sep = "")
  if (is.null(metrics) == T & NHmat_ans == F & genout_ans == 
      F) {
    stop("At least one metric must be specified or one of the following must be TRUE: NHmat_ans,\n         or genout_ans")
  }
  if (is.null(file_name) & NHmat_ans == T) {
    stop("If NHmat_ans = TRUE, you must specify a file_name")
  }
  if (is.null(file_name) & genout_ans == T) {
    stop("If genout_ans = TRUE, you must specify a file_name")
  }
  if (is.null(file_name) & "pFST" %in% metrics) {
    stop("If you include pFST as a metric, you must specify a file_name")
  }
  if (is.numeric(min_N) == F) {
    stop("min_N must be an integer")
  }
  if (is.numeric(xy[, 2]) == F) {
    stop("The second column of the xy file must be a numeric x coordinate (e.g. longitude)")
  }
  if (is.numeric(xy[, 3]) == F) {
    stop("The third column of the xy file must be a numeric y coordinate (e.g. latitude)")
  }
  if (nrow(xy) != N) {
    stop("The number of rows in the xy file does not match the number of individuals in the genepop file")
  }
  if (nrow(dist.mat) != N) {
    stop("The number of rows in the dist.mat does not match the number of individuals in the genepop file")
  }
  if (ncol(dist.mat) != N) {
    stop("The number of columns in the dist.mat does not match the number of individuals in the genepop file")
  }
  if (is.null(NeEstimator_dir) == F) {
    if (OS == "Windows") {
      if (file.exists(file.path(NeEstimator_dir, "Ne2.exe")) == 
          F) {
        stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")
      }
    }
    if (OS == "Linux") {
      if (file.exists(file.path(NeEstimator_dir, "Ne2L")) == 
          F) {
        stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")
      }
      if (file.exists(file.path(NeEstimator_dir, "NeEstimator.jar")) == 
          F) {
        stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")
      }
    }
    if (OS == "Darwin") {
      if (file.exists(file.path(NeEstimator_dir, "Ne2-1M")) == 
          F) {
        stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")
      }
      if (file.exists(file.path(NeEstimator_dir, "NeEstimator2x1.jar")) == 
          F) {
        stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")
      }
    }
  }
  if (is.null(NeEstimator_dir) == T & "NS" %in% metrics) {
    stop("You have specified NS as a metric, but you have not specified the location of NeEstimator_dir")
  }
  cat("Input summary:\n")
  cat(paste("\t individuals:", N, "\n"))
  cat(paste("\t loci:", numloci, "\n"))
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
  if (!is.null(file_name)) {
    sGD_outfile <- paste0(file_name, "_sGD.csv")
    cat(paste("\t sGD output file:", sGD_outfile, "in", getwd(), 
              "\n"))
  }
  if (NHmat_ans == TRUE) {
    NHmat_file <- paste0(file_name, "_neighborhood_matrix.csv")
    cat(paste("\t NH matrix file:", NHmat_file, "in", getwd(), 
              "\n"))
  }
  if (genout_ans == TRUE) {
    cat(paste("\t NH genepop file:", genout_file, "in", getwd(), 
              "\n"))
  }
  cat("Determining neighborhood membership from dist.mat and NH_radius...\n")
  NHmat <- c()
  for (r in c(1:N)) {
    if (length(NH_radius) == 1) {
      radius <- NH_radius
    }
    else {
      radius <- NH_radius[r]
    }
    rvals <- ifelse(dist.mat[r, ] < radius, 1, 0)
    NHmat <- rbind(NHmat, rvals)
  }
  NH_N <- as.numeric(rowSums(NHmat, na.rm = T))
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
  NH_summary <- data.frame(NH_Index, as.character(genind_obj@pop), 
                           xy[, c(2, 3)], NH_N, stringsAsFactors = F)
  names(NH_summary) <- c("NH_Index", "NH_ID", "X", "Y", "N")
  if (!is.null(metrics) | genout_ans == T) {
    header_length <- 1 + numloci
    genepop_header <- "sGD neighboorhood genepop file (each POP is a genetic neighborhood)"
    write.table(genepop_header, genout_file, sep = "\t", 
                quote = F, col.names = F, row.names = F)
    for (locus in 1:numloci) {
      write.table(names(genind_obj@all.names)[locus], genout_file, 
                  sep = "\t", quote = F, col.names = F, row.names = F, 
                  append = T)
    }
    for (pop in valid_pops) {
      write.table("POP", genout_file, sep = "\t", quote = F, 
                  col.names = F, row.names = F, append = T)
      if (!is.null(max_N) && NH_summary$N[pop] > max_N) {
        pop_ind <- sample(which(NHmat[pop, ] == 1), max_N)
      }
      else {
        pop_ind <- which(NHmat[pop, ] == 1)
      }
      NH_IDs <- paste0("POP", pop, "_", NH_summary$NH_ID[pop_ind])
      output <- df_dat[pop_ind, ]
      output <- data.frame(NH_IDs, ",", output[2:ncol(output)])
      #output[is.na(df_dat)] <- paste(rep("0", allele.digits * 2), collapse = "")
write.table(output, genout_file, sep = "\t", quote = F, 
            col.names = F, row.names = F, append = T)
    }
  }
  if ("GD" %in% metrics) {
    cat("Calculating gen,getic diversity indices...\n")
    NH_genind <- read.genepop(genout_file, ncode = allele.digits, 
                              quiet = T)
    #import2genind(genout_file, ncode=allele.digits)
    #NH_genind <- output
    NH_hierfstat <- genind2hierfstat(NH_genind)
    NH_A <- nb.alleles(NH_hierfstat)
    NH_Ar <- round(colMeans(allelic.richness(NH_hierfstat)$Ar), 
                   4)
    NH_Ap <- round(colSums(NH_A)/sum(NH_genind@loc.n.all), 
                   4)
    NH_A <- round(colMeans(NH_A), 4)
    NH_stats <- basic.stats(NH_hierfstat, digits = 4)
    NH_Ho <- round(colMeans(NH_stats$Ho), 4)
    NH_Hs <- round(colMeans(NH_stats$Hs), 4)
    NH_FIS <- round(colMeans(NH_stats$Fis), 4)
    GD_output <- data.frame(NH_Index[valid_pops], NH_A, NH_Ap, 
                            NH_Ar, NH_Hs, NH_Ho, NH_FIS, stringsAsFactors = F)
    names(GD_output) <- c("NH_Index", "A", "Ap", "Ar", "Hs", 
                          "Ho", "FIS")
    NH_summary <- merge(NH_summary, GD_output, by = "NH_Index", 
                        all = T)
  }
  if ("HWE" %in% metrics) {
    cat("Testing Hardy-Weinberg equilibrium for neighborhoods...\n")
    NH_diveRsity <- basicStats(genout_file)
    NH_HWE_p <- as.numeric(round(colMeans(NH_diveRsity$hwe_llr_p[, 
                                                                 -1], na.rm = T), 4))
    NH_HetEx <- as.numeric(round(colMeans(NH_diveRsity$hwe_het[, 
                                                               -1], na.rm = T), 4))
    NH_HomEx <- as.numeric(round(colMeans(NH_diveRsity$hwe_hom[, 
                                                               -1], na.rm = T), 4))
    HWE_output <- data.frame(NH_Index[valid_pops], NH_HWE_p, 
                             NH_HetEx, NH_HomEx, stringsAsFactors = F)
    names(HWE_output) <- c("NH_Index", "HWE_p", "HetEx_p", 
                           "HomEx_p")
    NH_summary <- merge(NH_summary, HWE_output, by = "NH_Index", 
                        all = T)
  }
  if ("NS" %in% metrics) {
    cat("Calculating NS...\n")
    write.table(1, "Ne2_input.txt", sep = "\t", col.names = F, 
                row.names = F, quote = F)
    write.table(3, "Ne2_input.txt", sep = "\t", col.names = F, 
                row.names = F, quote = F, append = T)
    write.table(cbind(0.1, 0.05, 0.02), "Ne2_input.txt", 
                sep = "\t", col.names = F, row.names = F, quote = F, 
                append = T)
    write.table(genout_file, "Ne2_input.txt", sep = "\t", 
                col.names = F, row.names = F, quote = F, append = T)
    write.table("Ne2_output.txt", "Ne2_input.txt", sep = "\t", 
                col.names = F, row.names = F, quote = F, append = T)
    if (OS == "Windows") {
      file.copy(file.path(NeEstimator_dir, "Ne2.exe"), 
                getwd())
      system("Ne2.exe  m:Ne2_input.txt", show.output.on.console = F, 
             ignore.stdout = T, ignore.stderr = T)
      file.remove("Ne2.exe")
    }
    if (OS == "Darwin") {
      file.copy(file.path(NeEstimator_dir, "Ne2-1M"), getwd())
      file.copy(file.path(NeEstimator_dir, "NeEstimator2x1.jar"), 
                getwd())
      system("./Ne2-1M  m:Ne2_input.txt", ignore.stdout = T, 
             ignore.stderr = T)
      file.remove("Ne2-1M")
      file.remove("NeEstimator2x1.jar")
    }
    if (OS == "Linux") {
      file.copy(file.path(NeEstimator_dir, "Ne2L"), getwd())
      file.copy(file.path(NeEstimator_dir, "NeEstimator 2.01.jar"), 
                getwd())
      system("./Ne2L  m:Ne2_input.txt", ignore.stdout = T, 
             ignore.stderr = T)
      file.remove("Ne2L")
      file.remove("NeEstimator 2.01.jar")
    }
    LDNe_output <- readLines("Ne2_output.txt")
    LDNe_datalines <- grep("Estimated Ne", LDNe_output)
    LDNe_data <- suppressWarnings(as.numeric(unlist(strsplit(LDNe_output[LDNe_datalines], 
                                                             "\\s+"))))
    LDNe_estimates <- data.frame(matrix(LDNe_data, nrow = npops, 
                                        ncol = 7, byrow = T)[, 4:7], stringsAsFactors = F)
    LDNe_estimates <- data.frame(NH_Index[valid_pops], LDNe_estimates, 
                                 stringsAsFactors = F)
    names(LDNe_estimates) <- c("NH_Index", "NS_ex0.10", "NS_ex0.05", 
                               "NS_ex0.02", "NS_ex0.00")
    NH_summary <- merge(NH_summary, LDNe_estimates, by = "NH_Index", 
                        all = T)
    file.remove("Ne2_input.txt")
    file.remove("Ne2_output.txt")
  }
  if (!is.null(file_name)) {
    cat("Writing sGD output file...\n")
    write.table(NH_summary, sGD_outfile, row.names = F, sep = ",", 
                na = "")
  }
  if (NHmat_ans == TRUE) {
    cat("Writing neighborhood membership matrix to file...\n")
    dimnames(NHmat) <- list(NH_summary$NH_ID, NH_summary$NH_ID)
    write.table(matrix(c(NA, NH_summary$NH_ID), nrow = 1), 
                NHmat_file, row.names = F, col.names = F, sep = ",", 
                na = "")
    write.table(NHmat, NHmat_file, row.names = NH_summary$NH_ID, 
                col.names = F, sep = ",", na = "", append = T)
  }
  if ("Dr" %in% metrics) {
    cat("Writing Roger's r genetic distance matrix to file...\n")
    if (!exists("NH_genind")) {
      NH_genind <- read.genepop(genout_file, ncode = allele.digits, 
                                quiet = T)
    }
    if (!exists("NH_hierfstat")) {
      NH_hierfstat <- genind2hierfstat(NH_genind)
    }
    Dr = genet.dist(NH_hierfstat, method = "Dr")
    write.table(full(Dr), paste0(file_name, "_Dr.csv"), sep = ",", 
                col.names = F, row.names = F)
  }
  if (genout_ans == FALSE) {
    file.remove(genout_file)
  }
  else {
    cat("Writing neighborhood genepop file...\n")
  }
  return(NH_summary)
  cat("Processing complete.\n")
  cat("\n")
}
  

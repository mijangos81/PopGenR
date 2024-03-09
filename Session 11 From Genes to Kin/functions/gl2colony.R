################################################################################
##          Function to export a COLONY2 file from a genlight object          ##
##                                                                            ##
##  Authors: Jesús Castrejón-Figueroa. R developer, Monash University         ##
##           Diana A Robledo-Ruiz. Research Fellow, Monash University         ##
##  Date: 2022-07-01                                                          ##
##                                                                            ##
##  This function requires:                                                   ##
##   - Input: a genlight object with offspring and parental information       ##
##            contained in gl@other$ind.metrics as columns 'offspring',       ##
##            'mother' and 'father' taking values 'yes'/'no'.                 ##
##   - User specified parameters:                                             ##
##         - filename_out = path of the output COLONY2 file.                  ##
##         - Others are set to a default value (check Colony User-Manual for  ##
##           further information).                                            ##
##   - Two custom made functions called 'parental.ids' and 'gl2structure'     ##
##     included after the main function in this script.                       ##
##                                                                            ##
##  Output:                                                                   ##
##   - A formatted text file ready to be analysed by COLONY2.                 ##
##                                                                            ##
##  Index:                                                                    ##
##    Line 31: Main function gl2colony                                        ##
##    Line 247: Function parental.ids                                         ##
##    Line 278: Function gl2structure                                         ##
##    Line 339: Example of use for gl2colony                                  ##
################################################################################

######################## Define MAIN function gl2colony ########################
library(crayon)

gl2colony <- function(gl = NULL,           
                      filename_out = NULL,
                      project_name = 'my_project',
                      output_name = 'my_project',
                      probability_father = 0.5,
                      probability_mother = 0.5,
                      seed = NULL,                # Seed for random number generator
                      update_allele_freq = 0,     # 0/1=Not updating/updating allele frequency
                      di_mono_ecious = 2,         # 2/1=Dioecious/Monoecious species
                      inbreed = 0,                # 0/1=no inbreeding/inbreeding
                      haplodiploid = 0,           # 0/1=Diploid species/HaploDiploid species
                      polygamy_male = 0,          # 0/1=Polygamy/Monogamy for males
                      polygamy_female = 0,        # 0/1=Polygamy/Monogamy for females
                      clone_inference = 1,        # 0/1=Clone inference =No/Yes
                      scale_shibship = 1,         # 0/1=Scale full sibship=No/Yes
                      sibship_prior = 0,          # 0/1/2/3/4=No/Weak/Medium/Strong/Optimal sibship prior
                      known_allele_freq = 0,      # 0/1=Unknown/Known population allele frequency
                      num_runs = 1,               # Number of runs
                      length_run = 2,             # 1/2/3/4=short/medium/long/very long run
                      monitor_method = 0,         # 0/1=Monitor method by Iterate#/Time in second
                      monitor_interval = 10000,   # Monitor interval in Iterate# / in seconds
                      windows_gui = 0,            # 0/1=No/Yes for run with Windows GUI
                      likelihood = 0,             # 0/1/2=PairLikelihood score/Fulllikelihood/FPLS
                      precision_fl = 2,           # 0/1/2/3=Low/Medium/High/Very high precision with Full-likelihood
                      marker_id = 'mk@',          # Marker Ids
                      marker_type = '0@',         # Marker types, 0/1=Codominant/Dominant
                      allelic_dropout = '0.000@', # Allelic dropout rate at each locus
                      other_typ_err = '0.05@',    # Other typing error rate at each locus
                      paternity_exclusion_threshold = '0 0',
                      maternity_exclusion_threshold = '0 0',
                      paternal_sibship = 0,
                      maternal_sibship = 0,
                      excluded_paternity = 0,
                      excluded_maternity = 0,
                      excluded_paternal_sibships = 0,
                      excluded_maternity_sibships = 0) {

  if(is.null(gl)){
      stop('Missing Genlight object.')
  }

  if(is.null(filename_out)){
      stop('Missing output filename.')
  }

  loci = dim(gl)[2]

  if(is.null(loci)){
      stop('Number of loci must be non-zero.')
  } 

  if(is.null(seed)){
    # random seed oherwise set manually
      seed <- sample.int(65535, 1)
  }
  message(sprintf('Random seed set to %d',seed))


  x <- parental.ids(gl)

  offspring_ids <- x$offs
  dad_ids       <- x$dad
  mum_ids       <- x$mum

  n_offspring <- length(offspring_ids)
  n_dads     <- length(dad_ids)
  n_mums   <- length(mum_ids)

  n_total <- n_offspring + n_mums + n_dads

  message(sprintf('%d Offspring detected. \n%d Fathers detected. \n%d Mothers detected.\n', n_offspring, n_dads, n_mums))

  if( !length(dad_ids) && !length(mum_ids)) {
    message('Missing paternal and maternal IDs: only sibship inference will be possible.')
  }

  if( !length(offspring_ids) ) {
    stop('Missing offspring IDs.')
  }

  # Transform genlight object to structure format
  message("Exporting Genlight object to COLONY2 format...")
  str_1row_with0s <- gl2structure(gl)

  # Subset structure dataset with offspring
  offspring_geno_to_keep <- match(offspring_ids, rownames(str_1row_with0s))  
  offspring_gen <- str_1row_with0s[offspring_geno_to_keep,]  # keep only offspring genotypes

  # Subset structure dataset with only mums
  if(length(mum_ids) > 0) {
    mum_geno_to_keep <- match(mum_ids, rownames(str_1row_with0s))
    mum_gen <- str_1row_with0s[mum_geno_to_keep,]  # keep only mum genotypes
  } else {
    probability_mother = 0
  }

  # Subset structure dataset with only dads
  if(length(dad_ids) > 0) {
    dad_geno_to_keep <- match(dad_ids, rownames(str_1row_with0s))
    dad_gen <- str_1row_with0s[dad_geno_to_keep,]  # keep only dad genotypes
  } else {
    probability_father = 0
  }

  ###################### 1. Create header for COLONY2 file
  head_comments <- list('! No. offspring',
  '! No. of loci',
  '! Seed for random number generator',
  '! 0/1=Not updating/updating allele frequency',
  '! 2/1=Dioecious/Monoecious species',
  '! 0/1=no inbreeding/inbreeding',
  '! 0/1=Diploid species/HaploDiploid species',
  '! 0/1=Polygamy/Monogamy for males & females',
  '! 0/1=Clone/duplicates inference =No/Yes',
  '! 0/1=Scale full sibship=No/Yes',
  '! 0/1/2/3=No/Weak/Medium/Strong sibship prior',
  '! 0/1=Unknown/Known population allele frequency',
  '! Number of runs',
  '! 1/2/3/4=short/medium/long/very long run',
  '! 0/1=Monitor method by Iterate#/Time in second',
  '! Monitor interval in Iterate# / in seconds',
  '! 0/1=No/Yes for run with Windows GUI',
  '! 0/1/2=PairLikelihood score/Fulllikelihood/FPLS',
  '! 0/1/2/3=Low/Medium/High/Very high precision FL',
  '',
  '! Marker Ids (consecutive for all)',
  '! Marker types, 0/1=Codominant/Dominant',
  '! Allelic dropout rate for all loci',
  '! Other typing error rate for all loci')

  polygamy <- paste(polygamy_male, polygamy_female,sep=' ')

  head.list <- list(n_offspring,loci,seed,update_allele_freq,di_mono_ecious,inbreed,haplodiploid,
                    polygamy,clone_inference,scale_shibship,sibship_prior,known_allele_freq,
                    num_runs,length_run,monitor_method,monitor_interval,windows_gui,likelihood,precision_fl,'',marker_id,
                    marker_type,allelic_dropout,other_typ_err
                    )
                    
  # Export header to COLONY2 file
  sink(filename_out)
  cat(project_name, '\n')
  cat(output_name, '\n')
  for(i in 1:length(head.list)) {
    cat(head.list[[i]], '\t\t', head_comments[[i]], '\n')
  }
  sink()
  
  ################# 2. Add offspring genotypes to COLONY2 file
  write.table(offspring_gen,
              file = filename_out,
              append = TRUE,
              quote = FALSE, 
              col.names = FALSE)

  sink(filename_out, append =TRUE)
  
  ################ 3. Add parents probabilities to COLONY2 file
  probabilities <- paste(probability_father, probability_mother, sep=' ')
  n_indv <- paste(n_dads, n_mums, sep=' ')
  cat('\n')
  cat(probabilities, '\t\t', 
      '! Probabilities that the father and mother of an offspring are included in candidates', 
      '\n')
  cat(n_indv, '\t\t', '! Numbers of candidate males and females', '\n')
  cat('\n')
  sink()
  
  ################# 4. Add dads genotypes to COLONY2 file
  if(length(dad_ids) > 0){
    message(sprintf("(%d%%) Working on it...", round(n_offspring*100/n_total,0)))
    write.table(dad_gen,
                file = filename_out,
                append = TRUE,
                quote = FALSE,
                col.names = FALSE)
  
    sink(filename_out, append = TRUE)
    cat('\n')
    sink()
  }
  
  ################# 5. Add mums genotypes to COLONY2 file
  if(length(mum_ids) > 0){
    message(sprintf("(%d%%) Almost there...", round((n_offspring + n_dads)*100/n_total,0)))
    write.table(mum_gen,
                file = filename_out,
                append = TRUE,
                quote = FALSE,
                col.names = FALSE)
  }
              
  ################ 6. Add last parameters to COLONY2 file
  last_comments <- list('! Number of offspring with known paternity, exclusion threshold',
                        '! Number of offspring with known maternity, exclusion threshold',
                        '',
                        '! Number of known paternal sibship',
                        '! Number of known maternal sibship',
                        '',
                        '! Number of offspring with known excluded paternity',
                        '! Number of offspring with known excluded maternity',
                        '',
                        '! Number of offspring with known excluded paternal sibships',
                        '! Number of offspring with known excluded maternal sibships')

  last_values <- list(paternity_exclusion_threshold,
                      maternity_exclusion_threshold,
                      '',
                      paternal_sibship,
                      maternal_sibship,
                      '',
                      excluded_paternity,
                      excluded_maternity,
                      '',
                      excluded_paternal_sibships,
                      excluded_maternity_sibships)

  sink(filename_out, append =TRUE)  
  cat('\n')
  
  for(i in 1:length(last_values)) {
    cat(last_values[[i]],'\t\t',last_comments[[i]],'\n')
  }
  sink()
  
  # Finished
  cat(crayon::green$bold('(100%) COLONY2 file successfully exported!'))
}
################################################################################


######################### Define function parental.ids #########################
## This function extracts parental information in a list of 3 elements (vectors 
## with offspring, dads and mums IDs, respectively)
parental.ids <- function(gen_data) {

  # Read metadata and convert to lowercase
  indv.metadata <- gen_data@other$ind.metrics
  names(indv.metadata) <- tolower(names(indv.metadata))

  # Remove leading/trailing white spaces
  indv.metadata$mother    <- tolower(indv.metadata$mother)
  indv.metadata$father    <- tolower(indv.metadata$father)
  indv.metadata$offspring <- tolower(indv.metadata$offspring)
  
  # Subset metadata
  mum_ids  <- indv.metadata[indv.metadata$mother    %in% c("yes"," yes", "yes "), 'id']
  dad_ids  <- indv.metadata[indv.metadata$father    %in% c("yes"," yes", "yes "), 'id']
  offs_ids <- indv.metadata[indv.metadata$offspring %in% c("yes"," yes", "yes "), 'id']
  
  # Make them vectors
  mum_ids  <- as.vector(na.omit(mum_ids))
  dad_ids  <- as.vector(na.omit(dad_ids))
  offs_ids <- as.vector(na.omit(offs_ids))

  # Make a list with the 3 vectors
  x = list(offs = offs_ids, dad = dad_ids, mum = mum_ids)
  return(x)
}
################################################################################


######################### Define function gl2structure #########################
## This function converts gl matrix to Structure format and from 2-row-per-ind 
## to 1-row-per-ind
gl2structure <- function(x,
                         addtlColumns = NULL, 
                         ploidy = 2,
                         exportMarkerNames = FALSE) {    

    genmat <- as.matrix(x)
    indNames <- dimnames(genmat)[[1]]
    nInd <- dim(genmat)[1] # number of individuals
    
    # Make sets of possible genotypes
    G <- list()
    for(i in 0:ploidy) {
        G[[i + 1]] <- c(rep(1, ploidy - i), rep(2, i))
    }
    #G[[ploidy + 2]] <- rep(-9, ploidy) # for missing data
    G[[ploidy + 2]] <- rep(0, ploidy) # for missing data

    # Set up data frame for Structure
    StructTab <- data.frame(ind = rep(indNames, each = ploidy))
    
    # Add any additional columns
    if(!is.null(addtlColumns)){
        for(i in 1:dim(addtlColumns)[2]){
            StructTab <- data.frame(StructTab, rep(addtlColumns[,i], each = ploidy))
            if(!is.null(dimnames(addtlColumns)[[2]])){
                names(StructTab)[i + 1] <- dimnames(addtlColumns)[[2]][i]
            } else {
                names(StructTab)[i + 1] <- paste("X", i, sep = "")
            }
        }
    }

    # Add genetic data
    for(i in 1:dim(genmat)[2]){
        thesegen <- genmat[, i] + 1
        thesegen[is.na(thesegen)] <- ploidy + 2
        StructTab[[dimnames(genmat)[[2]][i]]] <- unlist(G[thesegen])
    }
    
    # return(StructTab)  # Returning the value of gl2struct dartR function 
    
    data <- StructTab
    # Define dimensions of the matrix (only genotypes, not Ids)
    out <- matrix(NA, nrow = (nrow(data) / 2),        # no. of rows divided by 2
                      ncol = (2 * (ncol(data) - 1)))  # no. of columns minus Ids column times 2
  
    # Select first row per ind, leaving behind first column (Ids), then assign as first column per ind
    out[, seq(1, ncol(out), by = 2)] <- as.matrix(data[seq(1, nrow(data), by = 2), -1])
    # Select second row per ind, leaving behind first column (Ids), then assign as second column per ind
    out[, seq(2, ncol(out), by = 2)] <- as.matrix(data[seq(2, nrow(data), by = 2), -1])
  
    # Select Id column (only first row per ind) and make it rownames for matrix
    rownames(out) <- data[seq(1, nrow(data), by = 2), 1]  
    return(out) 
}
################################################################################

################################ Example of use ################################
##  gl2colony(gl = my.genlight,                                               ##
##            filename_out = "colony2.dat",                                   ##  
##            project_name = "parentage_fish_2022",                           ##  
##            output_name =  "parentage_fish_jul_2022",                       ##
##            seed = 1234,                                                    ##
##            probability_father = 0.6,                                       ##
##            probability_mother = 0.4,                                       ##
##            update_allele_freq = 1,                                         ##
##            allelic_dropout = '0.01@',                                      ##
##            other_typ_err = '0.001@')                                       ##
################################################################################

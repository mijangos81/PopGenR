
# function to sort genlight or genind alleles major/minor
# can specify pop(s) to use as reference if you want to analyse relative to some focal pop (optional)
sort_alleles <- function(geno, pops=NULL) {
  
  # Check if gl or gi supplied
  input.class <- class(geno)[1]
  
  # Convert genlight to genind if necessary
  wasGenlight <- FALSE
  if(input.class == "genlight") {
    geno <- gl2gi(geno)
    wasGenlight <- TRUE
  }
  
  
  # if focal population not supplied, calculate minor allele across all data
  if(is.null(pops)) {
    pops <- levels(geno@pop)
  }
  
  tab <- as.data.frame(geno@tab)
  tab$pop <- geno@pop
  tab <- tab[tab$pop %in% pops,]
  tab <- within(tab, rm(pop))
  
  res <- as.data.frame(colSums(tab, na.rm = T))
  colnames(res) <- "allele_count"
  res$allele_name <- rownames(res)
  
  
  my.ord <- do.call(c, lapply(seq(2, nrow(res), by = 2), (function(i){
    c(i-1, i)[order(res$allele_count[(i-1):i], decreasing = T)]})))
  
  
  res1 <- data.frame(res[my.ord,])
  
  geno@tab <- geno@tab[, res1$allele_name]
  
  # Convert back to genlight if the original input was genlight
  if(wasGenlight) {
    geno <- gi2gl(geno) # Assuming a hypothetical gi2gl function for conversion
  }
  
  return(geno)
  
}

# backward selection VIF function
vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  require(fmsb)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]))
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}


# R function rdadapt (see github.com/Capblancq/RDA-genome-scan)
rdadapt<-function(rda,K) {
  loadings<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(loadings, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

# function to expand pop-level metadata to individual-level dataframe
expand_pop2ind <- function(gl, meta.data) {
  # Extract individual names from genlight object
  ind_names <- indNames(gl)
  
  # Extract population assignments for each individual
  pops <- pop(gl)
  
  # Create a data frame with individual names and their population
  ind_df <- data.frame(ind.names = ind_names, site = pops, stringsAsFactors = FALSE)
  meta.data$site <- as.character(meta.data$site)
  
  # Use left_join to merge
  expanded_df <- ind_df %>%
    left_join(meta.data, by = "site")
  
  write.csv(expanded_df, "individual_env_data.csv", quote = FALSE, row.names = FALSE)
  return(expanded_df)
  
}


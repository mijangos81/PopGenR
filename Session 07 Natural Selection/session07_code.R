## ----warning=FALSE, message=FALSE-----------------------------------------------------------
#devtools::install_github("pygmyperch/melfuR")
#BiocManager::install('qvalue')
library(adegenet)
library(LEA)
library(vegan)
library(fmsb)
library(psych)
library(dartRverse)
library(melfuR)
library(sdmpredictors)
library(sf)
library(raster)
library(robust)
library(qvalue)
library(mlbench)
library(caret)
library(RAINBOWR)
library(qqman)
library(randomForest)
knitr::opts_knit$set(root.dir = "/cloud/project/")



## -------------------------------------------------------------------------------------------
source("/cloud/project/Session 07 Natural Selection/utils.R")


## -------------------------------------------------------------------------------------------
# load genlight object
load("/cloud/project/data/Mf5000_gl.RData")
Mf5000_gl
Mf5000_gl@pop

## What do you notice about the factor levels?
# Factors in R can make you question your life choices like nothing else...
# re-order the pop levels to match the order of individuals in your data
Mf5000_gl@pop <- factor(Mf5000_gl@pop, levels = as.character(unique(Mf5000_gl@pop)))
Mf5000_gl@pop

# convert to genind
Mf5000.genind <- gl2gi(Mf5000_gl)
Mf5000.genind

# Ok we have our genotypes


## ----echo=TRUE, eval=FALSE------------------------------------------------------------------
## # ADMIXTURE results most likely 6 pops
## imputed.gi <- melfuR::impute.data(Mf5000.genind, K = 6 )


## ----echo=FALSE, warning=FALSE, message=FALSE-----------------------------------------------
# suppress the 100s of lines of console output, instead print first 15 [...] last 5
capture_output <- capture.output({
imputed.gi <- melfuR::impute.data(Mf5000.genind, K = 6 )
})

output_length <- length(capture_output)

if (output_length > 20) {
cat(capture_output[1:15], sep="\n")
cat("\n\n[...]\n\n")
cat(capture_output[(output_length - 4):output_length], sep="\n")
} else {
cat(capture_output, sep="\n")
}


## -------------------------------------------------------------------------------------------
# re-order each genotype as major/minor allele
# (generally no real reason to do this... but it's neat and might be useful for something)
imputed.sorted.gi <- sort_alleles(imputed.gi)

# check results
imputed.gi@tab[1:10,1:6]
imputed.sorted.gi@tab[1:10,1:6]

# The allele counts for SNP_1 were ordered C/T before sorting, now they are T/C
# the sort_alleles function also has an optional 'pops' argument where you can specify
# a population to use as the reference if you want to analyse relative to some focal pop

unlink(c('dat.geno', 'dat.lfmm', 'dat.snmfProject', 'dat.snmf', 'dat.lfmm_imputed.lfmm',
         'individual_env_data.csv'), recursive = T)


## -------------------------------------------------------------------------------------------

# for individual based analyses
# get allele counts
alleles <- imputed.sorted.gi@tab

# get genotypes (counts of reference allele) and clean up locus names
snps <- alleles[,seq(1,ncol(alleles),2)]
colnames(snps) <- locNames(imputed.sorted.gi)
snps[1:10,1:10]


# Alternative: impute missing data with most common genotype across all data
# alleles <- Mf5000.genind@tab
# snps <- alleles[,seq(1,ncol(alleles),2)]
# colnames(snps) <- locNames(Mf5000.genind)
# snps <- apply(snps, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
# check total % missing data
# (sum(is.na(snps)))/(dim(snps)[1]*dim(snps)[2])*100
 


## -------------------------------------------------------------------------------------------
# for population based analyses
# get pop allele frequencies
gp <- genind2genpop(imputed.sorted.gi)
AF <- makefreq(gp)

# drop one (redundant) allele per locus
AF <- AF[,seq(1,ncol(AF),2)]
colnames(AF) <- locNames(gp)
rownames(AF) <- levels(imputed.sorted.gi@pop)
AF[1:14,1:6]



## ----echo=TRUE, eval=FALSE------------------------------------------------------------------
## # set the data directory and extend the response timeout (sometimes the database can be slow to respond)
## options(sdmpredictors_datadir="/cloud/project/data/spatial_data")
## options(timeout = max(300, getOption("timeout")))
## 
## # Explore environmental layers
## sdmpredictors::list_layers("WorldClim")
## 


## ----echo=FALSE, warning=FALSE, message=FALSE-----------------------------------------------
options(sdmpredictors_datadir="/cloud/project/data/spatial_data")
options(timeout = max(300, getOption("timeout")))

capture_output <- capture.output({
sdmpredictors::list_layers("WorldClim")
})

output_length <- length(capture_output)

if (output_length > 20) {
cat(capture_output[1:20], sep="\n")
cat("\n\n[...]\n\n")
} else {
cat(capture_output, sep="\n")
}


## -------------------------------------------------------------------------------------------
# Explore environmental layers
layers <- as.data.frame(sdmpredictors::list_layers("WorldClim"))
write.csv(layers, "/cloud/project/data/WorldClim.csv")

# Download specific layers to the datadir
# Bio01 (annual mean temperature), Bio05 (temperature in hottest month), Bio15 (rainfall seasonality) and Bio19 (rainfall in coldest quarter)
WC <- load_layers(c("WC_bio1", "WC_bio5", "WC_bio15", "WC_bio19"))
plot(WC)

# Crop rasters to study area and align CRS
# get Murray-Darling Basin polygon
mdb <- st_read("/cloud/project/data/spatial_data/MDB_polygon/MDB.shp")

# set CRS
mdb <- st_transform(mdb, 4326)

# convert to sf format
mdb <- as(mdb, "Spatial")

# crop WorldClim rasters to MDB extent, then mask pixels outside of the MDB
ENV <- crop(WC, extent(mdb))
ENV <- mask(ENV, mdb)


# get sampling sites, format as sf
sites <- read.csv("/cloud/project/data/Mf_xy.csv", header = TRUE)
sites_sf <- st_as_sf(sites, coords = c("X", "Y"), crs = 4326)


# Generate a nice color ramp and plot the rasters 
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

# include other spatial data that you might want to show
rivers <- st_read("/cloud/project/data/spatial_data/Mf_streams/Mf_streams.shp")
rivers <- st_transform(rivers, 4326)
rivers <- as(rivers, "Spatial")

# plot a layer
plot(ENV$WC_bio1, col=my.colors(100), main=names(ENV$WC_bio1))
lines(rivers, col="blue", lwd=0.3)
text(st_coordinates(sites_sf)[,1], st_coordinates(sites_sf)[,2], labels=sites_sf$site, cex=0.5)



## ----output=FALSE---------------------------------------------------------------------------

# format nicely and export as pdf
{pdf("/cloud/project/Session 07 Natural Selection/env_rasters.pdf")
for(i in 1:nlayers(ENV)) {
  rasterLayer <- ENV[[i]]
  
  # skew the colour ramp to improve visualisation (not necessary, but nice for data like this)
  p95 <- quantile(values(rasterLayer), probs = 0.05, na.rm = TRUE)
  breaks <- c(minValue(rasterLayer), seq(p95, maxValue(rasterLayer), length.out = 50))
  breaks <- unique(c(seq(minValue(rasterLayer), p95, length.out = 10), breaks))
  color.palette <- my.colors(length(breaks) - 1)
  layerColors <- if(i < 4) color.palette else rev(color.palette)
  
  # tidy up the values displayed in the legend
  simplifiedBreaks <- c(min(breaks), quantile(breaks, probs = c(0.25, 0.5, 0.75)), max(breaks))
  
  # plot the rasters
  plot(rasterLayer, breaks=breaks, col=layerColors, main=names(rasterLayer),
       axis.args=list(at=simplifiedBreaks, labels=as.character(round(simplifiedBreaks))),
       legend.args=list(text='', side=3, line=2, at=simplifiedBreaks))
  lines(rivers, col="blue", lwd=0.3, labels(rivers$Name))
  text(st_coordinates(sites_sf)[,1], st_coordinates(sites_sf)[,2], labels=sites_sf$site, cex=0.5)
}
dev.off()}



## -------------------------------------------------------------------------------------------

# and finally extract the environmental data for your sample sites
env.site <- data.frame(site = sites_sf$site,
                       X = st_coordinates(sites_sf)[,1],
                       Y = st_coordinates(sites_sf)[,2],
                       omegaX = sites_sf$omegaX,
                       omegaY = sites_sf$omegaY)

env.dat <- as.data.frame(raster::extract(ENV,st_coordinates(sites_sf)))
env.site <- cbind(env.site, env.dat)




## -------------------------------------------------------------------------------------------

# let's also just do a sanity check to make sure the genomic and environmental data are in the same order
stopifnot(all(env.site$site == rownames(AF)))



## -------------------------------------------------------------------------------------------
# to get a feel for the data we can run quick PCA using the rda function in vegan
pc <- rda(AF)
plot(pc)


## -------------------------------------------------------------------------------------------
# set colours, markers and labels for plotting and add to env.site df

env.site$cols <- c('darkturquoise', 'darkturquoise', 'darkturquoise', 'limegreen', 'limegreen', 'darkorange1', 'slategrey', 'slategrey', 'slategrey',
                   'slategrey', 'slategrey', 'dodgerblue3', 'gold', 'gold')

env.site$pch <- c(15, 16, 17, 17, 15, 15, 16, 15, 18, 16, 16, 8, 16, 16)

# generate labels for each axis
x.lab <- paste0("PC1 (", paste(round((pc$CA$eig[1]/pc$tot.chi*100),2)),"%)")
y.lab <- paste0("PC2 (", paste(round((pc$CA$eig[2]/pc$tot.chi*100),2)),"%)")



## ----eval=FALSE-----------------------------------------------------------------------------
## # plot PCA
## 
## {pdf(file = "/cloud/project/Session 07 Natural Selection/Mf_PCA.pdf", height = 6, width = 6)
## pcaplot <- plot(pc, choices = c(1, 2), type = "n", xlab=x.lab, ylab=y.lab, cex.lab=1)
## with(env.site, points(pc, display = "sites", col = env.site$cols, pch = env.site$pch, cex=1.5, bg = env.site$cols))
## legend("topleft", legend = env.site$site, col=env.site$cols, pch=env.site$pch, pt.cex=1, cex=0.50, xpd=1, box.lty = 0, bg= "transparent")
## dev.off()
## }


## ----include = FALSE------------------------------------------------------------------------
{png(file = "/cloud/project/Session 07 Natural Selection/Mf_PCA.png", units = 'in', height = 6, width = 6, res = 200)
  pcaplot <- plot(pc, choices = c(1, 2), type = "n", xlab=x.lab, ylab=y.lab, cex.lab=1)
  with(env.site, points(pc, display = "sites", col = env.site$cols, pch = env.site$pch, cex=1.5, bg = env.site$cols))
  legend("topleft", legend = env.site$site, col=env.site$cols, pch=env.site$pch, pt.cex=1, cex=0.50, xpd=1, box.lty = 0, bg= "transparent")
  dev.off()
}



## -------------------------------------------------------------------------------------------

# extract predictor variables to matrix
env_var <- env.site[ , c("WC_bio1", "WC_bio5", "WC_bio15", "WC_bio19")]

# check for obvious correlations between variables
pairs.panels(env_var, scale=T, lm = TRUE)

## reduce variance associated with covarying variables using variance inflation factor analyses
# run procedure to remove variables with the highest VIF one at a time until all remaining variables are below your threshold
keep.env <-vif_func(in_frame=env_var,thresh=2,trace=T)
keep.env  # the retained environmental variables
reduced.env <- subset(as.data.frame(env_var), select=c(keep.env))
scaled.env <- as.data.frame(scale(reduced.env))

# lets have another look
pairs.panels(reduced.env, scale=T, lm = TRUE)



## -------------------------------------------------------------------------------------------
# format spatial coordinates as matrix
xy <- as.matrix(cbind(env.site$omegaX, env.site$omegaY))



# reduced RDA model using the retained environmental PCs 
# conditioned on the retained spatial variables
Mf.RDA <- rda(AF ~ WC_bio5 + WC_bio15 + WC_bio19 + Condition(xy), data = scaled.env)
Mf.RDA


# So how much genetic variation can be explained by our environmental model?
RsquareAdj(Mf.RDA)$r.squared

# how much inertia is associated with each axis
screeplot(Mf.RDA)

# calculate significance of the reduced model

mod_perm <- anova.cca(Mf.RDA, nperm=1000) #test significance of the model
mod_perm

#generate x and y labels (% constrained variation)
x.lab <- paste0("RDA1 (", paste(round((Mf.RDA$CCA$eig[1]/Mf.RDA$CCA$tot.chi*100),2)),"%)")
y.lab <- paste0("RDA2 (", paste(round((Mf.RDA$CCA$eig[2]/Mf.RDA$CCA$tot.chi*100),2)),"%)")




## ----eval=FALSE-----------------------------------------------------------------------------
## #plot RDA1, RDA2
## {pdf(file = "/cloud/project/Session 07 Natural Selection/Mf_RDA.pdf", height = 6, width = 6)
## pRDAplot <- plot(Mf.RDA, choices = c(1, 2), type="n", cex.lab=1, xlab=x.lab, ylab=y.lab)
## with(env.site, points(Mf.RDA, display = "sites", col = env.site$cols, pch = env.site$pch, cex=1.5, bg = env.site$cols))
## text(Mf.RDA, "bp",choices = c(1, 2), labels = c("WC_bio5", "WC_bio15", "WC_bio19"), col="blue", cex=0.6)
## legend("topleft", legend = env.site$site, col=env.site$cols, pch=env.site$pch, pt.cex=1, cex=0.50, xpd=1, box.lty = 0, bg= "transparent")
## dev.off()}
## 


## ----include=FALSE--------------------------------------------------------------------------
#plot RDA1, RDA2
{png(file = "/cloud/project/Session 07 Natural Selection/Mf_RDA.png", units = 'in', res = 200, height = 6, width = 6)
pRDAplot <- plot(Mf.RDA, choices = c(1, 2), type="n", cex.lab=1, xlab=x.lab, ylab=y.lab)
with(env.site, points(Mf.RDA, display = "sites", col = env.site$cols, pch = env.site$pch, cex=1.5, bg = env.site$cols))
text(Mf.RDA, "bp",choices = c(1, 2), labels = c("WC_bio5", "WC_bio15", "WC_bio19"), col="blue", cex=0.6)
legend("topleft", legend = env.site$site, col=env.site$cols, pch=env.site$pch, pt.cex=1, cex=0.50, xpd=1, box.lty = 0, bg= "transparent")
dev.off()}



## -------------------------------------------------------------------------------------------

k = 2 # number of RDA axes to retain
rda.q <- rdadapt(Mf.RDA, K = k)
sum(rda.q$q.values <= 0.01)
candidate.loci <- colnames(AF)[which(rda.q[,2] < 0.01)]


# check the distribution of candidates and q-values
ggplot() +
  geom_point(aes(x=c(1:length(rda.q[,2])), 
                 y=-log10(rda.q[,2])), col="gray83") +
  geom_point(aes(x=c(1:length(rda.q[,2]))[which(rda.q[,2] < 0.01)], 
                 y=-log10(rda.q[which(rda.q[,2] < 0.01),2])), col="orange") +
  xlab("SNPs") + ylab("-log10(q.values") +
  theme_bw()


# check correlations between candidate loci and environmental predictors
# this section relies on code by Brenna Forester, found at
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html

npred <- 3
env_mat <- matrix(nrow=(sum(rda.q$q.values < 0.01)), 
                  ncol=npred)  # n columns for n predictors
colnames(env_mat) <- colnames(reduced.env)


# calculate correlation between candidate snps and environmental variables
for (i in 1:length(candidate.loci)) {
  nam <- candidate.loci[i]
  snp.gen <- AF[,nam]
  env_mat[i,] <- apply(reduced.env,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(as.data.frame(candidate.loci),env_mat)  

for (i in 1:length(cand$candidate.loci)) {
  bar <- cand[i,]
  cand[i,npred+2] <- names(which.max(abs(bar[,2:4]))) # gives the variable
  cand[i,npred+3] <- max(abs(bar[,2:4]))              # gives the correlation
}

colnames(cand)[npred+2] <- "predictor"
colnames(cand)[npred+3] <- "correlation"

table(cand$predictor) 
write.csv(cand, "/cloud/project/data/RDA_candidates.csv", row.names = FALSE)



## -------------------------------------------------------------------------------------------
env.ind <- expand_pop2ind(Mf5000_gl, env.site)

# now subset our genotype objects to the candidate loci
candidate.gl <- gl.keep.loc(gi2gl(imputed.sorted.gi), loc.list=candidate.loci)
candidate.gi <- gl2gi(candidate.gl)

# get genotypes (counts of reference allele) and clean up locus names
candidate.alleles <- candidate.gi@tab
candidate.snps <- candidate.alleles[,seq(1,ncol(candidate.alleles),2)]
colnames(candidate.snps) <- locNames(candidate.gi)

# another quick sanity check
stopifnot(all(env.ind$ind.names == indNames(candidate.gi)))


# finally run the individual based RDA
# no need to control for pop structure this time
Mfcandidate.RDA <- rda(candidate.snps ~ WC_bio5 + WC_bio15 + WC_bio19, data = env.ind)
Mfcandidate.RDA

ind.mod_perm <- anova.cca(Mfcandidate.RDA, nperm=1000, 
                          parallel=4) #test significance of the model
ind.mod_perm

x.lab <- paste0("RDA1 (",
                paste(round((Mfcandidate.RDA$CA$eig[1]/Mfcandidate.RDA$tot.chi*100),2)),"%)")
y.lab <- paste0("RDA2 (",
                paste(round((Mfcandidate.RDA$CA$eig[2]/Mfcandidate.RDA$tot.chi*100),2)),"%)")




## ----eval=FALSE-----------------------------------------------------------------------------
## 
## #plot RDA1, RDA2
## {pdf(file = "/cloud/project/Session 07 Natural Selection/Mfcandidate_RDA.pdf", height = 6, width = 6)
## pRDAplot <- plot(Mfcandidate.RDA, choices = c(1, 2),
##                  type="n", cex.lab=1, xlab=x.lab, ylab=y.lab)
## points(Mfcandidate.RDA, display = "sites", col = env.ind$cols,
##        pch = env.ind$pch, cex=0.8, bg = env.ind$cols)
## text(Mfcandidate.RDA, "bp",choices = c(1, 2),
##      labels = c("WC_bio5", "WC_bio15", "WC_bio19"), col="blue", cex=0.6)
## legend("topleft", legend = env.site$site, col=env.site$cols, pch=env.site$pch,
##        pt.cex=1, cex=0.50, xpd=1, box.lty = 0, bg= "transparent")
## dev.off()
## }
## 
## 


## ----include = FALSE------------------------------------------------------------------------

#plot RDA1, RDA2
{png(file = "/cloud/project/Session 07 Natural Selection/Mfcandidate_RDA.png", units = 'in', res = 200, height = 6, width = 6)
pRDAplot <- plot(Mfcandidate.RDA, choices = c(1, 2), type="n", cex.lab=1, xlab=x.lab, ylab=y.lab)
points(Mfcandidate.RDA, display = "sites", col = env.ind$cols, pch = env.ind$pch, cex=0.8, bg = env.ind$cols)
text(Mfcandidate.RDA, "bp",choices = c(1, 2), labels = c("WC_bio5", "WC_bio15", "WC_bio19"), col="blue", cex=0.6)
legend("topleft", legend = env.site$site, col=env.site$cols, pch=env.site$pch, pt.cex=1, cex=0.50, xpd=1, box.lty = 0, bg= "transparent")
dev.off()
}





## -------------------------------------------------------------------------------------------

## Get input files 
load(file = "/cloud/project/data/GWAS.RData")
##randomForest
CaRF$Size <- as.factor(CaRF$Size)



## -------------------------------------------------------------------------------------------

See(CaRF)

table(CaRF$Size)



## -------------------------------------------------------------------------------------------

set.seed(123)
ind <- sample(2, nrow(CaRF), replace = TRUE, prob = c(0.75, 0.25))
train <- CaRF[ind==1,]
test <- CaRF[ind==2,]


rf <- randomForest(Size~., data=train,
                   ntree = 50,       # number of trees
                   mtry = 10,          # number of variables tried at each split
                   importance = TRUE,
                   proximity = TRUE)
print(rf)



## -------------------------------------------------------------------------------------------

control <- trainControl(method="oob", number=5)
set.seed(123) 
metric <- "Accuracy"
mtry <- 10
tunegrid <- expand.grid(.mtry=mtry )
rf_small <- train(Size~., data=CaRF, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control, ntree=50)
rf_small



## -------------------------------------------------------------------------------------------

control <- trainControl(method="oob", number=5)
set.seed(123) 
metric <- "Accuracy"
mtry <- sqrt(ncol(CaRF))
tunegrid <- expand.grid(.mtry=mtry )
rf_default <- train(Size~., data=CaRF, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
rf_default



## ----eval = FALSE---------------------------------------------------------------------------
## control <- trainControl(method="oob", number=10, search="random")
## set.seed(123)
## mtry <- sqrt(ncol(CaRF))
## rf_random <- train(Size~., data=CaRF, method="rf", metric=metric, tuneLength=15, trControl=control)
## (rf_random)
## plot(rf_random)


## -------------------------------------------------------------------------------------------
####Expected output
###To Change 
#  mtry  Accuracy   Kappa    
#195  0.7703453  0.5403670
#348  0.7685185  0.5367943
#526  0.7684685  0.5367281
#665  0.7721972  0.5441283
#953  0.7730731  0.5458833
#1011  0.7694444  0.5386467
#1017  0.7702202  0.5401206
#1038  0.7685435  0.5367490
#1115  0.7740240  0.5477867
#1142  0.7750000  0.5497308
#1253  0.7740490  0.5478857
#1268  0.7740240  0.5477339
#1627  0.7768018  0.5533787
#1842  0.7768268  0.5533826
#2013  0.7749750  0.5496527
####


## ----eval = FALSE---------------------------------------------------------------------------
## 
## control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
## tunegrid <- expand.grid(.mtry= c(565, 665, 765))
## modellist <- list()
## seed<-123
## metric <- "Accuracy"
## for (ntree in c(200, 300, 500, 1000)) {
##   set.seed(seed)
##   fit <- train(Size~., data=CaRF, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control, ntree=ntree)
##   key <- toString(ntree)
##   modellist[[key]] <- fit
## }
## # Then you can compare tuning results
## results <- resamples(modellist)
## summary(results)
## dotplot(results)


## -------------------------------------------------------------------------------------------

control <- trainControl(method="oob", number=5)
set.seed(123) 
metric <- "Accuracy"
mtry <- 665
tunegrid <- expand.grid(.mtry=mtry)
rf_best <- train(Size~., data=CaRF, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control, ntree=200)
rf_best



## -------------------------------------------------------------------------------------------

varImpPlot(rf_best$finalModel)
varImpPlot(rf)



## -------------------------------------------------------------------------------------------

modify.CaGW.res <- modify.data(pheno.mat = CaGW_phe, geno.mat = CaGW_geno, map = CaGW_map, return.ZETA = TRUE, return.GWAS.format = TRUE)
pheno.GWAS <- modify.CaGW.res$pheno.GWAS
geno.GWAS <- modify.CaGW.res$geno.GWAS

### View each dataset
See(pheno.GWAS)  
See(geno.GWAS)
See(CaGW_map)



## -------------------------------------------------------------------------------------------

ZETA <- modify.CaGW.res$ZETA



## ----eval = FALSE---------------------------------------------------------------------------
## 
## ### Population structure, using a structure Q values matrix, n.PC should be > 0
## normal.res <- RGWAS.normal(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA, structure.matrix= QMATRIX, n.PC = 4, P3D = TRUE)
## 
## #### other cofactors, like sex or age (a two-column table similar to phenotype table)
## normal.res <- RGWAS.normal(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA, structure.matrix= QMATRIX, covariate = SEX_AGE,  n.PC = 4, P3D = TRUE)


## ----results='hide', message=FALSE----------------------------------------------------------

normal.res <- RGWAS.normal(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA, n.PC = 4, P3D = TRUE)




## -------------------------------------------------------------------------------------------

normal.res$D[normal.res$D$Size  >= 4, ]



## -------------------------------------------------------------------------------------------

K.A2<-calcGRM(genoMat = CaGW_geno) ### aditive effects
k.AA2 <- K.A2 * K.A2 ### epistatic effects

y <- as.matrix(pheno.GWAS[,'Size', drop = FALSE])
Z <- design.Z(pheno.labels = rownames(y), geno.names = rownames(K.A2))
pheno.mat2 <- y[rownames(Z), , drop = FALSE]
ZETA2 <- list(A = list(Z = Z, K = K.A2),AA = list(Z = Z, K = k.AA2))
EM3.res <- EM3.cpp(y = pheno.mat2, X = NULL, ZETA = ZETA2) ### multi-kernel linear mixed effects model
(Vu <- EM3.res$Vu)   ### estimated genetic variance
(Ve <- EM3.res$Ve)   ### estimated residual variance
(weights <- EM3.res$weights) ### estimated proportion of two genetic variances
(herit <- Vu * weights / (Vu + Ve)) ### genomic heritability
herit



## ----eval = FALSE---------------------------------------------------------------------------
## 
## phenoNoNA <- pheno.mat2[noNA, , drop = FALSE]   ### remove NA
## 


## -------------------------------------------------------------------------------------------

nFold <- 10 
nLine <- nrow(pheno.mat2)
idCV <- sample(1:nLine %% nFold) ### assign random ids for cross-validation
idCV[idCV == 0] <- nFold

yPred <- rep(NA, nLine)
for (noCV in 1:nFold) {
  print(paste0("Fold: ", noCV))
  yTrain <- pheno.mat2
  yTrain[idCV == noCV, ] <- NA ### prepare test data
  EM3.gaston.resCV <- EM3.general(y = yTrain, X0 = NULL, ZETA = ZETA,
                                  package = "gaston", return.u.always = TRUE,
                                  pred = TRUE, return.u.each = TRUE,
                                  return.Hinv = TRUE) 
  yTest <- EM3.gaston.resCV$y.pred ### predicted values
  yPred[idCV == noCV] <- yTest[idCV == noCV]
}



## -------------------------------------------------------------------------------------------

plotRange <- range(pheno.mat2, yPred)
plot(x = pheno.mat2, y = yPred,xlim = plotRange, ylim = plotRange,
     xlab = "Observed values", ylab = "Predicted values",
     main = "Results of Genomic Prediction (multi-kernel)",
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.3)
abline(a = 0, b = 1, col = 2, lwd = 2, lty = 2)
R2 <- cor(x = pheno.mat2[, 1], y = yPred) ^ 2
text(x = plotRange[2] - 10,
     y = plotRange[1] + 10,
     paste0("R2 = ", round(R2, 3)),
     cex = 1.5)



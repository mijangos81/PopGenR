## ----warning=FALSE, message=FALSE--------------------------------------
library(dartRverse)
library(ggplot2)
knitr::opts_knit$set(root.dir = "/cloud/project")

## ----------------------------------------------------------------------

# load data
load("./data/Session3_data.RData")

#create a list of the kimberley genlights
kimberley_names <- ls(pattern = "^Kimberley")
#put all the genlights into a mega list
kimberley <- mget(kimberley_names)

#we're going to do this in a loop for speed, applying the same filters
# Iterate over the names of the kimberley list
for(name in names(kimberley)){
  # Extract the genlight object from the kimberley list using its name
  genlight_object <- kimberley[[name]]
  # Apply the filter call rate function
  # Assuming gl.filter.callrate is a function that operates on a genlight object
  filtered_object <- gl.filter.callrate(genlight_object, threshold = 0.7, mono.rm = TRUE)
  # Assign the filtered object back to the environment with a new name
  assign(paste0(name, "_0.7"), filtered_object)
}


## ----------------------------------------------------------------------
# List all object names in the environment
all_names <- ls()

# Use grep() to match names that start with "Kimberley" 
# and end with "0.7", .+ indicates any characters in between
kimberley_filtered <- grep("^Kimberley.+0\\.7$", all_names, value = TRUE)#put all the genlights into a mega list
#create another list with the ones we want
kimberley <- mget(kimberley_filtered)

#Initialize an empty data frame
heterozygosity_reports_df <- data.frame()

# Iterate over the kimberley list to apply gl.report.heterozygosity
# and bind the results
for(name in names(kimberley)) {
  # Apply the function
  report <- gl.report.heterozygosity(kimberley[[name]])
  
  # Add 'ObjectName' as the first column of the report
  report <- cbind(ObjectName = name, report)
  
  # Bind this report to the main data frame
  heterozygosity_reports_df <- bind_rows(heterozygosity_reports_df, report)
}

# heterozygosity_reports_df now contains all the reports with an 
# additional column for object names
heterozygosity_reports_df



## ----------------------------------------------------------------------
kimberley_Ho_0.7callrate <- ggplot(heterozygosity_reports_df, aes(x = ObjectName, y = Ho)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  labs(title = "Observed Heterozygosity by Sample number", x = "Sample Number", y = "Observed Heterozygosity (Ho) at 0.7 Call Rate")

kimberley_Ho_0.7callrate


## ----------------------------------------------------------------------
# reload data
load('./data/Session3_data.RData')
# List all object names in the environment
all_names <- ls()

# Use grep() to match names that start with "Kimberley" 
kimberley_names <- grep("^Kimberley.*\\.vcf$", all_names, value = TRUE) #put all the genlights into a mega list
#create another list with the ones we want
kimberley <- mget(kimberley_names)

#we're going to do this in a loop for speed, applying the same filters
# Iterate over the names of the kimberley list
for(name in names(kimberley)){
  # Extract the genlight object from the kimberley list using its name
  genlight_object <- kimberley[[name]]
  # Apply the filter call rate function
  # Assuming gl.filter.callrate is a function that operates on a genlight object
   filtered_object <- gl.filter.callrate(genlight_object, threshold = .95, mono.rm = TRUE)
  # Assign the filtered object back to the environment with a new name
  assign(paste0(name, "_0.95"), filtered_object)
}


## ----------------------------------------------------------------------
# List all object names in the environment
all_names <- ls()

# Use grep() to match names that start with "Kimberley" 
# and end with "0.95", .+ indicates any characters in between
kimberley_filtered <- grep("^Kimberley.+0\\.95$", all_names, value = TRUE)#put all the genlights into a mega list
#create another list with the ones we want
kimberley <- mget(kimberley_filtered)


#Initialize an empty data frame
heterozygosity_reports_df_0.95 <- data.frame()

# Iterate over the kimberley list to apply gl.report.heterozygosity
# and bind the results
for(name in names(kimberley)) {
  # Apply the function
  report <- gl.report.heterozygosity(kimberley[[name]])
  
  # Add 'ObjectName' as the first column of the report
  report <- cbind(ObjectName = name, report)
  
  # Bind this report to the main data frame
  heterozygosity_reports_df_0.95 <- bind_rows(heterozygosity_reports_df_0.95, report)
}

# heterozygosity_reports_df now contains all the reports with an 
# additional column for object names
heterozygosity_reports_df_0.95


## ----------------------------------------------------------------------

# Example using ggplot2 to plot the data
library(ggplot2)
kimberley_Ho_0.95callrate <- ggplot(heterozygosity_reports_df_0.95, aes(x = ObjectName, y = Ho)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  labs(title = "Observed Heterozygosity by Sample number 0.95 Call Rate",
       x = "Sample Number", y = "Observed Heterozygosity (Ho) at 0.95 Call Rate")
kimberley_Ho_0.95callrate

par(mfrow = c(2,1))
kimberley_Ho_0.7callrate +kimberley_Ho_0.95callrate



## ----------------------------------------------------------------------
# Use grep() to match names that start with "Kimberley" 
kimberley_names <- grep("^Kimberley.*\\.vcf$", all_names, value = TRUE) 
# put all the genlights into a mega list

#create another list with the ones we want
kimberley <- mget(kimberley_names)

# Now lets subsample the datasets down to the same five individuals
# so that the only difference is our SNP calling
# who are the individuals
inds <- indNames(Kimberley_n_05.vcf_0.7)

#Initialize an empty data frame
heterozygosity_when_subsampling <- data.frame()

for(name in names(kimberley)) {
  # Access the genlight object from your list
  genlight_object <- kimberley[[name]]
  
  # Subset the individuals
  x <- gl.keep.ind(genlight_object, ind.list = inds, mono.rm = TRUE)
  
  #filter on call rate
  
  x <- gl.filter.callrate(x, threshold = .7)
  
  # Apply the function
  report <- gl.report.heterozygosity(x)
  
  # Add 'ObjectName' as the first column of the report
  report$ObjectName <- name 
  # Adjusting to add column without cbind to maintain data frame classes
  
  # Bind this report to the main data frame
  heterozygosity_when_subsampling <- bind_rows(heterozygosity_when_subsampling,
                                               report)
}



## ----------------------------------------------------------------------
kimberley_Ho_subsampling <- ggplot(heterozygosity_when_subsampling, aes(x = ObjectName, y = Ho)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  labs(title = "Observed Heterozygosity when subsampling", 
       x = "SNP Calling", y = "Observed Heterozygosity (Ho)")

kimberley_Ho_subsampling



## ----------------------------------------------------------------------
#create a list of the southeast genlights
Southeast_names <- ls(pattern = "^SouthEast")
#put all the genlights into a mega list
southeast <- mget(Southeast_names)

#Now lets subsample the datasets down to the same five individuals
#so that the only difference is our SNP calling
#who are the individuals
inds <- indNames(SouthEast_n_05.vcf)


for(name in names(southeast)) {
  # Access the genlight object from your list
  genlight_object <- southeast[[name]]
  
  # Subset the individuals
  x <- gl.keep.ind(genlight_object, ind.list = inds, mono.rm = TRUE)
  
  #filter on call rate
  
  x <- gl.filter.callrate(x, threshold = .7)
  
  # Apply the function
  report <- gl.report.heterozygosity(x)
  
  # Add 'ObjectName' as the first column of the report
  report$ObjectName <- name # Adjusting to add column without cbind to maintain data frame classes
  
  # Bind this report to the main data frame
  heterozygosity_when_subsampling <- bind_rows(heterozygosity_when_subsampling, report)
}

southeast_Ho_subsampling <- ggplot(heterozygosity_when_subsampling, aes(x = ObjectName, y = Ho)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  labs(title = "Observed Heterozygosity when subsampling",
       x = "SNP Calling", y = "Observed Heterozygosity (Ho)")


## ----------------------------------------------------------------------
#create a list of the central australian genlights
central_names <- ls(pattern = "^Central")
#put all the genlights into a mega list
central <- mget(central_names)

#Now lets subsample the datasets down to the same five individuals
#so that the only difference is our SNP calling
#who are the individuals
inds <- indNames(Central_n_05.vcf)

for(name in names(central)) {
  # Access the genlight object from your list
  genlight_object <- central[[name]]
  
  # Subset the individuals
  x <- gl.keep.ind(genlight_object, ind.list = inds, mono.rm = TRUE)
  
  #filter on call rate
  
 x <- gl.filter.callrate(x, threshold = .7)
  
  # Apply the function
  report <- gl.report.heterozygosity(x)
  
  # Add 'ObjectName' as the first column of the report
  report$ObjectName <- name # Adjusting to add column without cbind to maintain data frame classes
  
  # Bind this report to the main data frame
  heterozygosity_when_subsampling <- bind_rows(heterozygosity_when_subsampling, report)
}

central_Ho_subsampling <- ggplot(heterozygosity_when_subsampling, aes(x = ObjectName, y = Ho)) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
  labs(title = "Observed Heterozygosity when subsampling 5 individuals",
       x = "SNP Calling", y = "Observed Heterozygosity (Ho)")




## ----------------------------------------------------------------------
southeast_Ho_subsampling
central_Ho_subsampling


## ----------------------------------------------------------------------
#first we assign populations to the individuals in the genlight
individual_names <- as.vector(indNames(combined_kim_15_cen_15_se_15.vcf))

# Split each name at the underscore and keep the first part
modified_names <- sapply(individual_names, function(name) {
  parts <- strsplit(name, "_")[[1]]
  parts[1]
})

# Convert the output to a factor
modified_names <- factor(modified_names)

#then we assign the populations
combined_kim_15_cen_15_se_15.vcf@pop <- modified_names

#then we apply a call rate filter
combined_kim_15_cen_15_se_15.vcf_0.7 <- gl.filter.callrate(combined_kim_15_cen_15_se_15.vcf,
                                                           threshold = 0.95, mono.rm = T)

#then we look at the heterozygosity estimates by population
diversity_equal_sample_15_ind <- gl.report.heterozygosity(combined_kim_15_cen_15_se_15.vcf_0.7)


## ----------------------------------------------------------------------
individual_names <- as.vector(indNames(combined_kim_15_cen_5_se_10.vcf))

# Split each name at the underscore and keep the first part
modified_names <- sapply(individual_names, function(name) {
  parts <- strsplit(name, "_")[[1]]
  parts[1]
})

# Convert the output to a factor
modified_names <- factor(modified_names)

#then we assign the populations
combined_kim_15_cen_5_se_10.vcf@pop <- modified_names

#then we apply a call rate filter
combined_kim_15_cen_5_se_10.vcf_0.95 <- gl.filter.callrate(combined_kim_15_cen_5_se_10.vcf, threshold = 0.95, mono.rm = T)

#then we look at the heterozygosity estimates by population
diversity_unequal_sample <- gl.report.heterozygosity(combined_kim_15_cen_5_se_10.vcf_0.95)

#how much more heterozygosity does kimberley have than central in different analyses?

#where our samples are equal?
diversity_equal_sample_15_ind$Ho[1]/diversity_equal_sample_15_ind$Ho[3]

#where our samples are equal?
diversity_unequal_sample$Ho[1]/diversity_unequal_sample$Ho[3]


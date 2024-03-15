
#installs the necessary bioconductor packages
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("SNPRelate")

#install dartRverse (dartRverse) & core (dartR.base, dartR.data)
install.packages("dartRverse")

#now install all the packages you want

dartRverse::install("dartR.popgen") #cran version
dartRverse::install("dartR.popgen", "github", "main") #github version, main branch (tested)

dartRverse::install("dartR.popgen", "github", "dev") #dev version


#to check which packages and versions are installed

dartRverse::install()


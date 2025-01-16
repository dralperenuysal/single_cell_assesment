library(GEOquery)

wd <- "/home/alperen/PycharmProjects/single_cell_assesment/"

# set the working directory
setwd(wd)

# define the GEO number
geo <- "GSE266577"

# get the gse information
gse <- getGEO(GEO = geo,
              GSEMatrix = T)
gse <- gse[[1]] # there is only one single platform

# get the phenotype data and store it for further analysis
phenotype <- pData(phenoData(gse))

# save the phenotype data
write.table(phenotype,
            file = paste0(wd, "resource/phenotype.tsv"),
            sep = '\t')

# supplementary files downloaded manually because the GEOquery packages failed.

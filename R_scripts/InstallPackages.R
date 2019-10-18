#This short R script installs R packages required for DESeq2 analysis.
#Packages: DESeq2, RColorBrewer, gplots, getopt, optparse, data.table

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("RColorBrewer")
BiocManager::install("gplots")

install.packages("https://cran.r-project.org/src/contrib/Archive/getopt/getopt_1.20.0.tar.gz", repos=NULL, type="source", lib=".")
install.packages("https://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.4.4.tar.gz", repos=NULL, type="source", lib=".")
install.packages("https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.10.4.tar.gz", repos=NULL, type="source", lib=".")

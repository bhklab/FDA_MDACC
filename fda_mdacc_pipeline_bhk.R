########################
## Benjamin Haibe-Kains
## All rights Reserved
## July 15, 2015
########################

## remove all existing objects from the workspace
# rm(list=ls(all=TRUE))

## choose your favorite mirror
chooseCRANmirror(graphics=FALSE, ind=15)

## set path to local directory if it is not properly set up
.libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))

## load and install libraries
library(parallel)
library(elasticnet)
library(pvclust)
library(Hmisc)
library(vcd)
library(epibasix)
library(plotrix)
library(gdata)
library(WriteXLS)
library(tools)

## install the last version of genefu and survcomp
library(devtools)
## survcomp
availib <- require(survcomp)
if (!availib) {
  devtools::install_github(repo="bhklab/survcomp", ref="master")
  library(survcomp)
}
## genefu
availib <- require(genefu)
if (!availib) {
  devtools::install_github(repo="bhklab/genefu", ref="master")
  library(genefu)
}

## define functions required for the analysis pipeline
source(file.path("fda_mdacc_foo_bhk.R"))

## directory where all the analysis results will be stored
saveres <- "output"
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE, recursive=TRUE) }


## path to the raw data
rawpath <- "/mnt/work1/users/bhklab/Data/FDA/MDACC2"
## note that the files have been compressed using gzip (.CEL.gz)

## file name for normalized dataset
dataset.fn <- file.path(saveres, "mdacc_2_data_frma.RData")

########################
## global parameters

## set method for downloading
# options(download.file.method="auto")
options(download.file.method="wget")
## change to curl, wget or internal depending on your system

## prevent strings to be converted into factors
options(stringsAsFactors=FALSE)

## set random seed to ensuer reproducibility of the resuls
set.seed(54321)

## number of cpu cores available for the analysis pipeline
## set to 'NULL' if all the available cores should be used
nbcore <- 8
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore }
options("mc.cores"=nbcore)

## list of characters to be removed from row and column names
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"



########################
## run analysis scripts

message("\n---------------------------------------\n| Compute cutoff for GGI and PI3KCASIG |\n---------------------------------------")
source(file.path("fda_mdacc_ggi_cutoff.R"))
source(file.path("fda_mdacc_pi3kcasig_cutoff.R"))
message("\t-> DONE")

message("\n-----------------------------------\n| Format and normalize dataset |\n-----------------------------------")
source(file.path("fda_mdacc_dataset_bhk.R"))
message("\t-> DONE")

message("\n-----------------------------------\n| Compute signature |\n-----------------------------------")
source(file.path("fda_mdacc_signature_bhk.R"))
message("\t-> DONE")

## session info
writeLines(capture.output(sessionInfo()), file.path(saveres, "sessionInfo.log"))


## end
  


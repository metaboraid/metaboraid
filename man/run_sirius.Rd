% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_sirius.R
\name{run_sirius}
\alias{run_sirius}
\title{Perform identification using FINGER-ID}
\usage{
run_sirius(
  parameter_zip_files = NA,
  database = NA,
  local_database = NA,
  number_of_compounds_for_ion = 1,
  number_of_compounds = 10,
  use_heuristic = T,
  mz_to_use_heuristic_only = 650,
  mz_to_use_heuristic = 300,
  ncores = 1,
  progress_bar = F,
  verbose = F,
  conda = "auto",
  env = "metaboraid_package",
  results_folder = NA,
  chech_file_interval = 2,
  total_check_time = 20,
  timeout = 600,
  timeout_tree = 0
)
}
\arguments{
\item{parameter_zip_files}{A zip file containing MS2 parameter generated by \code{\link{map_adducts}}}

\item{database}{Name of the database to use: One of all, chebi, kegg, bio, natural products, pubmed, hmdb, biocyc, hsdb, knapsack, biological, zinc bio, gnps, pubchem, mesh, and maconda!}

\item{number_of_compounds_for_ion}{Number of compounds to consider for building trees. Default 1}

\item{number_of_compounds}{Number of compounds to report for each ion. Default 10}

\item{use_heuristic}{The heuristic will be used to initially rank the formula candidates. Default TRUE}

\item{mz_to_use_heuristic_only}{Set minimum m/z to only use heuristic tree computation. Default 650}

\item{mz_to_use_heuristic}{Set minimum m/z to enable heuristic preprocessing. Default 300}

\item{ncores}{Number of cores to use for parallel processing. Default 1}

\item{progress_bar}{Whether to show progress bar or not. Default FALSE}

\item{verbose}{Show information about different stages of the processes. Default FALSE}

\item{conda}{Conda binary to use. Default auto}

\item{env}{conda environment used to run the process. Default metaboraid_package}

\item{results_folder}{a path to a folder where results for EACH ion will be returned}

\item{chech_file_interval}{not used}

\item{total_check_time}{not used}

\item{timeout}{The maximum number of seconds to wait for a single out to give a result. Default 600}

\item{timeout_tree}{the number of seconds to wait for tree computation 0= unlimited. Default 0}
}
\value{
A dataframe containing the identified ions. The dataframe contains search engine and database specific information but also tree important columns: parentMZ, parentRT, fileName which are used to trace the ions by the downstream processes.
}
\description{
This function allows you to perform identification using FINGER-ID/
Please note that you need to have FINGER-ID, e.g using \code{\link{install_tools}}
}
\details{
We install CSI-FINGERID in metaboraid_package environment.
At this stage, CSI-FINGERID, Metfrag and Chemdistiller can be installed on Linux, IOS and Windows however, you should install CFM-ID yourself (\url{https://cfmid.wishartlab.com/})
}
\examples{

library(CAMERA)
library(metaboraid)
# Read MS1 and MS2 files
ms1_files<-system.file("ms1data",c("X1_Rep1.mzML","X2_Rep1.mzML"),package = "metaboraid")
ms2_files<-system.file("ms2data",c("sample1.mzML","sample2.mzML"),package = "metaboraid")

# mass trace detection
xs <- xcmsSet(ms1_files,method="centWave",ppm=30,peakwidth=c(5,10))

# mass trace matching
xsg <- group(xs)

# convert to CAMERA
xsa <- xsAnnotate(xsg)

# Group mass traces
anF <- groupFWHM(xsa, perfwhm = 0.6)

# Detect isotopes
anI <- findIsotopes(anF, mzabs = 0.01)

# Group using correlation
anIC <- groupCorr(anI, cor_eic_th = 0.75)

# Find adducts
anFA <- findAdducts(anIC, polarity="positive")

# map features and MS2s
mapped_features<-map_features(inputMS2s = ms2_files,input_camera = anFA,ppm = 10,rt = 10)

# Map adducts
mapped_adducts<-map_adducts(inputMS2List=mapped_features,input_camera=anFA,
                            precursorppm=10,
                            fragmentppm=20,fragmentabs=0.01,minPrecursorMass=NA,maxPrecursorMass=NA,
                            minPeaks=10,maxSpectra=10,mode="pos",adductRules="primary",
                            outputDir="general_parameters_4",searchMultipleChargeAdducts=T,
                            includeMapped=T,includeUnmapped=F,verbose=T)

# Run the search

run_sirius("parameter_files.zip",database = "all",ncores = 10,progress_bar = T,verbose = T,results_folder = "pp",chech_file_interval = 2,timeout = 600,conda = "auto")

}

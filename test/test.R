
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
ID_results<-run_sirius("parameter_files.zip",database = "all",ncores = 2,progress_bar = T,verbose = T,results_folder = "results",chech_file_interval = 2,timeout = 600,conda = "auto")

# map ID to features

final_results<-metaboraid:::map_to_camera(input_camera = anFA,identification_result = ID_results,ppm = 10,rt = 10,higher_the_better = T,score_column = "SiriusScore",impute = T,mapping_method = "fast")

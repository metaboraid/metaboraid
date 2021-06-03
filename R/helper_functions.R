#' Create a random directory with prefix in R temp dir
#'
#' Especially for testing code it is very helpful to have
#' a temp directory with a defined prefix, so one knows which test
#' produced which directory.
#' @param prefix A string that will preceed the directory name
#' @param dir Directory where the random dir will be created under.
#'            Defaults to tempdir()
#' @return Name of temporary directory
#' @keywords internal
createTmpDir <- function(prefix=NULL, dir=tempdir()) {
  tmpname <- NULL
  while(is.null(tmpname)){
    tmpname <- try(
      {
        tmp <- file.path(dir, paste(c(prefix, sample(letters,8)), collapse=""))
        makeDir(tmp, overwrite="never")
        return(tmp)
      }, silent=TRUE)
    if (class(tmpname)=="try-error") {
      tmpname <- NULL # restarts while
    }
  }
  return(tmpname)
}


#' Make a directory after performing an existence check
#'
#' Throws an exception if file or directory with same name exist
#' and overwrite is TRUE.
#' @param dir Name of directory to create
#' @param overwrite A character string: never (default), erase, overwrite
#' @return Path to created directory
#' @keywords internal
makeDir <- function(dir, overwrite="never"){
  if (file.exists(dir)) {
    if (overwrite=="never") stop("io.R/makeDir: I won't overwrite data present in dir=", dir)
    if (overwrite=="erase") safeUnlink(dir)
  }
  dir.create(dir, recursive=TRUE, showWarnings=FALSE)
  return(dir)
}

#' Checks if a particular software is installed and throws an exception it does not
#' @param software_name Name of the software
#' @param env conda environment to run the software
#' @param conda conda binary
#' @return if failded it will return 1
#' @keywords internal
#' @import reticulate
check_if_software_exist<-function(software_name=NA,env="metaboraid_package",conda="auto")
{
  if(!is.na(software_name) & !is.na(software_name))
  {
    out <- tryCatch({system2(reticulate::conda_binary(conda),shQuote(c("run","-n",env,software_name)),stdout = T,stderr = T)},
                    error=function(cond) {
                      message(sprintf("Probably %s has not been installed properly!",software_name))
                      return(1)
                    }, error=function(cond) {})
    return(out)
  }
}

#' Checks if a package on conda has been installed
#' @param pkg Name of the package
#' @param conda conda binary
#' @return none
#' @keywords internal
#' @import reticulate
check_if_conda_is_installed<-function(pkg="metaboraid_package",conda="auto")
{

  if (is_windows()) {

    # conda is the only supported method on windows
    method <- "conda"

    # confirm we actually have conda
    have_conda <- !is.null(tryCatch(reticulate::conda_binary(conda), error = function(e) NULL))
    if (!have_conda) {
      stop("Software installation failed (no conda binary found)\n\n",
           "Install Anaconda for Python 3.x (https://www.anaconda.com/download/#windows)\n",
           "before installing",
           call. = FALSE)
    }

    # avoid DLL in use errors
    if (py_available()) {
      stop("You should call install_tools() only in a fresh ",
           "R session that has not yet initialized Keras and TensorFlow (this is ",
           "to avoid DLL in use errors during installation)")
    }
  }else{
    have_conda <- !is.null(tryCatch(reticulate::conda_binary(conda), error = function(e) NULL))
    if (!have_conda) {
      stop("Software installation failed (no conda binary found)\n\n",
           "Install Anaconda for Python 3.x (https://www.anaconda.com/download/)\n",
           "before running this function",
           call. = FALSE)
    }
  }
  list_of_env<-reticulate::conda_list(conda)
  if(!pkg%in%list_of_env[,"name"])
  {
    stop(sprintf("We did not find the %s environment. Use install_tools to install the tools first!",pkg))
  }
}

#' Check if Windows is the running OS
#' @return a logical value TRUE if Windows is the OS
#' @keywords internal
is_windows <- function() {
  identical(.Platform$OS.type, "windows")
}

#' Check if IOS is the running OS
#' @return a logical value TRUE if IOS is the OS
#' @keywords internal
is_osx <- function() {
  Sys.info()["sysname"] == "Darwin"
}


#' Correct strangely looking adducts
#'
#' @param adduct a string containing addcut
#' @return the corrected addcut if needed correction
#' @keywords internal
validate.adduct <- function(adduct) {
  if(adduct == "[M+H+NH3]+") {return("[M+NH4]+")}
  if(adduct == "M+Na-2H") {return("M-2H+Na")}
  if(adduct == "M+K-2H") {return("M-2H+K")}
  return (adduct)
}

#' Generates metfrag format parameters
#'
#' Given the list of parameters, this function will generate a file containing those paramters
#'
#' @param param Parameter list
#' @param outputName Name a file for output
#' @return none
#' @keywords internal
#'
parameterToCommand<-function(param,outputName="")
{
  param$MetFragPeakListReader<-"de.ipbhalle.metfraglib.peaklistreader.FilteredStringTandemMassPeakListReader"
  param$MetFragCandidateWriter<-"CSV"
  param$PeakListString<-paste(apply(param$PeakList,1,paste,collapse="_"),collapse = ";")
  param$SampleName<-outputName
  param[[which(names(param)=="PeakList")]]<-NULL

  param$MetFragDatabaseType
  toOutput<-""
  for(i in 1:length(param))
  {
    if(i==1) {toOutput<- paste(names(param)[i],"=",param[[i]],sep="")
    } else {toOutput<- paste(toOutput," ",names(param)[i],"=",param[[i]],sep="")}
  }
  ## the output is specid_rt_mz_intensity[_origfilename]. specid is an enumerative number and origfilename is read from the original mzml file if included
  cat(toOutput, file = outputName, sep="\n")
}

#' Generate addcut information for MS2 data
#'
#' Given the list of MS2s (mapped and/or unmapped), this function estimate the neutral mass and attempts to associate addcuts to MS2 ions
#'
#' @param mappedMS2 A mapped MS2 list generated by \code{\link{map_features}}
#' @param unmappedMS2 A unmapped MS2 list generated by \code{\link{map_features}}
#' @param cameraObject A CAMERA object (results of preprocessing using CAMERA)
#' @param minPrecursorMass Minimum precursor mass to consider. Default 0
#' @param maxPrecursorMass Maximum precursor mass to consider NA=all. Default NA
#' @param minPeaks The minimum number of MS2 peaks that should be present for an ion to be consider
#' @param maxSpectra The maximum number of ions to write (mostly used for testing)
#' @param mode Ionization mode. Either pos or neg corresponding to positive and negative, respectively
#' @param primary What types of adducts to consider. If true primary otherwise extended (see details). Default primary
#' @param searchMultipleChargeAdducts Set to TRUE if you want to consider the selected adducts (primary or extended) for unmapped ions or ions without adduct information from CAMERA. Default TRUE
#' @param includeMapped Set to TRUE if you want to include mapped ions. Default TRUE
#' @param includeUnmapped Set to TRUE if you want to include unmapped ions. Default FALSE
#' @param settingsObject a list. Can be empty or build using another round of mapping on another file
#' @return a list of parameters
#' @details
#' MetaboRAID attempts to estimate the neutral mass of the molecules to decrease the search space of search engines.
#' The approach starts with mapping the MS2 parent mz and RT to the CAMERA features. We then search whether that specific feature has been annotated with adducts or isoptopes.
#' If so, the neutral mass will be estimated based on the adduct information.
#' However if there is no adduct information in CAMERA results for a particular ion, if searchMultipleChargeAdducts is set to TRUE, we generate all the selected ions based on adductRules:
#' primary: [M-H]-, [M-2H+Na]-, [M-2H+K]-, [M+Cl]-, [M+H]+, [M+Na]+, [M+K]+, [M+NH4]+
#' extended: We genrated up to three changes state for adduct information for each of them
#'
#' @keywords internal
#'
#' @import CAMERA
#' @import stringr
#' @import MSnbase
toMetfragCommand<-function(mappedMS2=NA,
                           unmappedMS2=NA,
                           cameraObject=NA,
                           searchMultipleChargeAdducts=F,
                           includeUnmapped=T,includeMapped=T,
                           settingsObject=list(),preprocess=NA,savePath="",minPeaks=0,maxSpectra=NA,
                           maxPrecursorMass = NA, minPrecursorMass = NA, mode = "pos", primary = T)
{
  peakList<-CAMERA::getPeaklist(cameraObject)
  file.origin<-""
  numberSpectraWritten <- 0
  if(includeMapped)
  {
    searchChargeFlag<-F
    searchAdductsFlag<-F
    massTraceNames<-names(mappedMS2)
    metFragResult<-c()
    for(x in massTraceNames)
    {
      seachAdducts<-NA
      seachCharge<-NA
      searchChargeFlag<-F
      intb <- peakList[as.numeric(x),"intb"]
      if(peakList[as.numeric(x),"adduct"]=="" & peakList[as.numeric(x),"isotopes"]=="")
      {
        adduct<-NA
        neutralMASS<-peakList[as.numeric(x),"mz"]
        searchChargeFlag<-T
        searchAdductsFlag<-T
        seachAdducts<-adduct
      }else if(peakList[as.numeric(x),"adduct"]!="")
      {
        if(stringr::str_count(peakList[as.numeric(x),"adduct"], "]")==1)
        {
          adduct<-gsub(" ","",stringr::str_extract(peakList[as.numeric(x),"adduct"], "\\[.*\\]*. "))
          neutralMASS<-as.numeric(stringr::str_extract(peakList[as.numeric(x),"adduct"], " .*"))
          searchChargeFlag<-F
        }else
        {
          adduct<-""
          neutralMASS<-peakList[as.numeric(x),"mz"]
          adduct<-unique(sapply(strsplit(stringr::str_extract_all(peakList[as.numeric(x),"adduct"], "\\[.*?\\]")[[1]]," "),
                                function(x){x[[1]]}))
          searchChargeFlag<-T
          seachAdducts<-adduct
        }
      }else if(peakList[as.numeric(x),"isotopes"]!="")
      {

        isotopID<- stringr::str_extract(peakList[as.numeric(x),"isotopes"], "\\[\\d*\\]")

        monoIsotopic<-grepl("[M]",peakList[ grepl(isotopID,peakList[,"isotopes"],fixed=T),"isotopes"],fixed=T)
        adduct<- gsub(" ","",
                      stringr::str_extract(peakList[ grepl(isotopID,peakList[,"isotopes"],fixed=T),][monoIsotopic,"adduct"],"\\[.*\\]*. "))

        tmpMASS<-as.numeric(stringr::str_extract(peakList[ grepl(isotopID,peakList[,"isotopes"],fixed=T),][monoIsotopic,"adduct"], " .*"))
        if(adduct!="" & !is.na(adduct)[1])
        {
          if(stringr::str_count(peakList[ grepl(isotopID,peakList[,"isotopes"],fixed=T),][monoIsotopic,"adduct"], "]")==1)
          {
            adduct<-gsub(" ","",stringr::str_extract(peakList[ grepl(isotopID,peakList[,"isotopes"],fixed=T),][monoIsotopic,"adduct"], "\\[.*\\]*. "))
            neutralMASS<-as.numeric(stringr::str_extract(peakList[ grepl(isotopID,peakList[,"isotopes"],fixed=T),][monoIsotopic,"adduct"], " .*"))
            searchChargeFlag<-F
          }else
          {
            adduct<-""
            neutralMASS<-peakList[as.numeric(x),"mz"]
            adduct<-unique(sapply(strsplit(stringr::str_extract_all(peakList[ grepl(isotopID,peakList[,"isotopes"],fixed=T),][monoIsotopic,"adduct"], "\\[.*?\\]")[[1]]," "),
                                  function(x){x[[1]]}))
            searchChargeFlag<-T
            seachAdducts<-adduct
          }

        }else
        {
          isoTMP<-gsub("\\[.*\\]","",peakList[as.numeric(x),"isotopes"])
          charges<-stringr::str_extract(isoTMP,"\\d+")

          if(is.na(charges))
          {
            adduct<-NA
            neutralMASS<-peakList[as.numeric(x),"mz"]
            searchChargeFlag<-T
            seachCharge<-1
          }else
          {
            neutralMASS<-peakList[as.numeric(x),"mz"]
            searchChargeFlag<-T
            seachCharge<-charges
          }

        }

      }

      mappedMS2TMP<-NA
      if(class(mappedMS2[[x]])=="list")
      {
        mappedMS2TMP<-mappedMS2[[x]]
      } else if(class(mappedMS2[[x]])=="Spectrum2")
      {
        mappedMS2TMP<-list(mappedMS2[[x]])
      }
      for(MSMS in mappedMS2TMP)
      {
        if(preprocess)
        {

          MSMS@centroided<-F
          MSMS@polarity<-as.integer(1)
          MSMS@smoothed<-F
          MSMS<-MSnbase::pickPeaks(MSMS)

        }
        MS2<-as.matrix(cbind(MSMS@mz,MSMS@intensity))
        # if number MS/MS peaks is too low
        if(length(MSMS@mz) == 0) { next }
        if(!is.na(minPeaks) & dim(MS2)[1] < minPeaks) { next }
        if(!searchChargeFlag)
        {
          settingsObject[["NeutralPrecursorMass"]]<-neutralMASS
          settingsObject[["PeakList"]]<-MS2
          settingsObject[["IsPositiveIonMode"]]<-"True"
          if(mode == "neg") {settingsObject[["IsPositiveIonMode"]]<-"False"}
          modeSuffix<-"+"
          if(mode == "neg") {modeSuffix<-"-"}
          settingsObject[["PrecursorIonType"]]<-validate.adduct(adduct)
          fileName<-""
          # add id, rt, neu_mass, intensity, orig file name
          file.origin <- gsub("\\?.*", "", gsub(".*/", "", attributes(MSMS)$fileName))
          if(file.origin == "") {
            fileName<-paste(as.character(numberSpectraWritten+1),"_",as.character(MSMS@rt),"_",as.character(round(MSMS@precursorMz,4)),"_",as.character(intb),".txt",sep="")
          } else {
            fileName<-paste(as.character(numberSpectraWritten+1),"_",as.character(MSMS@rt),"_",as.character(round(MSMS@precursorMz,4)),"_",as.character(intb),"_",file.origin,".txt",sep="")
          }
          if(savePath!="")
            fileName<-paste(savePath,"/",fileName,sep="")
          if(!is.na(maxPrecursorMass) & maxPrecursorMass < neutralMASS) { next }
          if(!is.na(minPrecursorMass) & minPrecursorMass > neutralMASS) { next }
          if(is.na(maxSpectra) || maxSpectra > numberSpectraWritten) {
            parameterToCommand(settingsObject,fileName)
            numberSpectraWritten<-numberSpectraWritten+1
          }
        } else if(searchChargeFlag & searchMultipleChargeAdducts)
        {

          allChargesHits<-list()
          allAdductForSearch<-adductCalculator(mz = neutralMASS,charge = seachCharge,
                                               adduct = gsub("\\[|\\]","",seachAdducts),mode = mode,primary = primary)
          for(k in 1:nrow(allAdductForSearch))
          {
            mass <- allAdductForSearch[k,"correctedMS"]
            settingsObject[["NeutralPrecursorMass"]]<-mass
            settingsObject[["PeakList"]]<-MS2
            settingsObject[["IsPositiveIonMode"]]<-"True"
            if(mode == "neg") {settingsObject[["IsPositiveIonMode"]]<-"False"}
            modeSuffix<-"+"
            if(mode == "neg") {modeSuffix<-"-"}
            settingsObject[["PrecursorIonType"]]<-paste("[",validate.adduct(as.character(allAdductForSearch[k,"adductName"])),"]", modeSuffix, sep="")
            fileName<-""
            file.origin <- gsub("\\?.*", "", gsub(".*/", "", attributes(MSMS)$fileName))
            if(file.origin == "") {
              fileName<-paste(as.character(numberSpectraWritten+1),"_",as.character(MSMS@rt),"_",as.character(round(MSMS@precursorMz,4)),"_",as.character(intb),".txt",sep="")
            } else {
              fileName<-paste(as.character(numberSpectraWritten+1),"_",as.character(MSMS@rt),"_",as.character(round(MSMS@precursorMz,4)),"_",as.character(intb),"_",file.origin,".txt",sep="")
            }
            if(savePath!="")
              fileName<-paste(savePath,"/",fileName,sep="")
            if(!is.na(maxPrecursorMass) & maxPrecursorMass < mass) { next }
            if(!is.na(minPrecursorMass) & minPrecursorMass > mass) { next }
            if(is.na(maxSpectra) || maxSpectra > numberSpectraWritten) {
              parameterToCommand(settingsObject,fileName)
              numberSpectraWritten<-numberSpectraWritten+1
            }
          }
        }
      }
    }
  }

  if(includeUnmapped)
  {
    for(p in 1:length(unmappedMS2))
    {
      MSMS<-unmappedMS2[[p]]
      if(preprocess)
      {

        MSMS@centroided<-F
        MSMS@polarity<-as.integer(1)
        MSMS@smoothed<-F
        MSMS<-MSnbase::pickPeaks(MSMS)

      }
      neutralMASS<-MSMS@precursorMz
      MS2<-as.matrix(cbind(MSMS@mz,MSMS@intensity))
      adduct<-"[M+H]+"
      if(mode == "neg") {adduct<-"[M-H]-"}
      if(length(MSMS@mz) == 0) { next }
      if(!is.na(minPeaks) & dim(MS2)[1] < minPeaks) { next }
      if(!searchMultipleChargeAdducts)
      {
        settingsObject[["NeutralPrecursorMass"]]<-neutralMASS
        settingsObject[["PeakList"]]<-MS2
        settingsObject[["IsPositiveIonMode"]]<-"True"
        if(mode == "neg") {settingsObject[["IsPositiveIonMode"]]<-"False"}
        settingsObject[["PrecursorIonType"]]<-adduct
        fileName<-""
        intb<-MSMS@precursorIntensity
        if(file.origin == "") {
          fileName<-paste(as.character(numberSpectraWritten+1),"_",as.character(MSMS@rt),"_",as.character(round(MSMS@precursorMz,4)),"_",as.character(intb),".txt",sep="")
        } else {
          fileName<-paste(as.character(numberSpectraWritten+1),"_",as.character(MSMS@rt),"_",as.character(round(MSMS@precursorMz,4)),"_",as.character(intb),"_",file.origin,".txt",sep="")
        }
        if(savePath!="")
          fileName<-paste(savePath,"/",fileName,sep="")
        if(!is.na(maxPrecursorMass) & maxPrecursorMass < neutralMASS) { next }
        if(!is.na(minPrecursorMass) & minPrecursorMass > neutralMASS) { next }
        if(is.na(maxSpectra) || maxSpectra > numberSpectraWritten) {
          parameterToCommand(settingsObject,fileName)
          numberSpectraWritten<-numberSpectraWritten+1
        }
      }else if(searchMultipleChargeAdducts)
      {
        allChargesHits<-list()
        allAdductForSearch<-adductCalculator(mz = neutralMASS,charge = NA,
                                             adduct = NA,mode = mode, primary = primary)
        for(k in 1:nrow(allAdductForSearch))
        {
          mass <- allAdductForSearch[k,"correctedMS"]
          settingsObject[["NeutralPrecursorMass"]]<-mass
          settingsObject[["PeakList"]]<-MS2
          settingsObject[["IsPositiveIonMode"]]<-"True"
          if(mode == "neg") {settingsObject[["IsPositiveIonMode"]]<-"False"}
          modeSuffix<-"+"
          if(mode == "neg") {modeSuffix<-"-"}
          settingsObject[["PrecursorIonType"]]<-paste("[",as.character(allAdductForSearch[k,"adductName"]),"]", modeSuffix, sep="")
          fileName<-""
          if(file.origin == "") {
            fileName<-paste(as.character(numberSpectraWritten+1),"_",as.character(MSMS@rt),"_",as.character(round(MSMS@precursorMz,4)),"_",as.character(intb),".txt",sep="")
          } else {
            fileName<-paste(as.character(numberSpectraWritten+1),"_",as.character(MSMS@rt),"_",as.character(round(MSMS@precursorMz,4)),"_",as.character(intb),"_",file.origin,".txt",sep="")
          }
          if(savePath!="")
            fileName<-paste(savePath,"/",fileName,sep="")
          if(!is.na(maxPrecursorMass) & maxPrecursorMass < mass) { next }
          if(!is.na(minPrecursorMass) & minPrecursorMass > mass) { next }
          if(is.na(maxSpectra) || maxSpectra > numberSpectraWritten) {
            parameterToCommand(settingsObject,fileName)
            numberSpectraWritten<-numberSpectraWritten+1
          }
        }
      }
    }

  }
}


#' Generate neutral mass for different addcuts
#'
#' Given mz, charge, ionization, adduct information this fucntion will estimate neutral mass for addcuts
#'
#' @param mz A mz value
#' @param charge possible charge of the ion
#' @param mode ionization mode. Can be pos or neg
#' @param adduct possible addcut
#' @param primary If TRUE, will use primary adducts otherwise up to 3 chanrges will be used
#' @details
#' if searchMultipleChargeAdducts is set to TRUE, we generate all the selected ions based on adductRules:
#' primary: [M-H]-, [M-2H+Na]-, [M-2H+K]-, [M+Cl]-, [M+H]+, [M+Na]+, [M+K]+, [M+NH4]+
#' extended: We genrated up to three changes state for adduct information for each of them
#'
adductCalculator<-function(mz=NA,charge=NA,mode="pos",adduct=NA, primary = T)
{
  if(is.na(mz))stop("Provide the mz!")
  if(is.na(mode))stop("Provide the polarity!")

  Ion.name.pos=c("M+3H",
                 "M+2H+Na",
                 "M+H+2Na",
                 "M+3Na",
                 "M+2H",
                 "M+H+NH4",
                 "M+H+Na",
                 "M+H+K",
                 "M+ACN+2H",
                 "M+2Na",
                 "M+2ACN+2H",
                 "M+3ACN+2H",
                 "M+H",
                 "M+NH4",
                 "M+Na",
                 "M+CH3OH+H",
                 "M+K",
                 "M+ACN+H",
                 "M+2Na-H",
                 "M+IsoProp+H",
                 "M+ACN+Na",
                 "M+2K-H",
                 "M+DMSO+H",
                 "M+2ACN+H",
                 "M+IsoProp+Na+H",
                 "2M+H",
                 "2M+NH4",
                 "2M+Na",
                 "2M+K",
                 "2M+ACN+H",
                 "2M+ACN+Na")

  Ion.mass.pos=c(
    "M/3 + 1.007276",
    "M/3 + 8.334590",
    "M/3 + 15.7661904",
    "M/3 + 22.989218",
    "M/2 + 1.007276",
    "M/2 + 9.520550",
    "M/2 + 11.998247",
    "M/2 + 19.985217",
    "M/2 + 21.520550",
    "M/2 + 22.989218",
    "M/2 + 42.033823",
    "M/2 + 62.547097",
    "M + 1.007276",
    "M + 18.033823",
    "M + 22.989218",
    "M + 33.033489",
    "M + 38.963158",
    "M + 42.033823",
    "M + 44.971160",
    "M + 61.06534",
    "M + 64.015765",
    "M + 76.919040",
    "M + 79.02122",
    "M + 83.060370",
    "M + 84.05511",
    "2M + 1.007276",
    "2M + 18.033823",
    "2M + 22.989218",
    "2M + 38.963158",
    "2M + 42.033823",
    "2M + 64.015765")

  Charge.pos=c(
    "3+",
    "3+",
    "3+",
    "3+",
    "2+",
    "2+",
    "2+",
    "2+",
    "2+",
    "2+",
    "2+",
    "2+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+",
    "1+")

  Extended.pos<-c(T,T,T,T,T,T,T,T,T,T,T,T,F,F,F,T,F,T,T,T,T,T,T,T,T,T,T,T,T,T,T)

  Ion.name.neg<-c("M-3H",
                  "M-2H",
                  "M-H2O-H",
                  "M-H",
                  "M+Na-2H",
                  "M+Cl",
                  "M+K-2H",
                  "M+FA-H",
                  "M+Hac-H",
                  "M+Br",
                  "M+TFA-H",
                  "2M-H",
                  "2M+FA-H",
                  "2M+Hac-H",
                  "3M-H")

  Ion.mass.neg<-c("M/3 - 1.007276",
                  "M/2 - 1.007276",
                  "M - 19.01839",
                  "M - 1.007276",
                  "M + 20.974666",
                  "M + 34.969402",
                  "M + 36.948606",
                  "M + 44.998201",
                  "M + 59.013851",
                  "M + 78.918885",
                  "M + 112.985586",
                  "2M - 1.007276",
                  "2M + 44.998201",
                  "2M + 59.013851",
                  "3M - 1.007276")

  Charge.neg<-c("3-",
                "2-",
                "1-",
                "1-",
                "1-",
                "1-",
                "1-",
                "1-",
                "1-",
                "1-",
                "1-",
                "1-",
                "1-",
                "1-",
                "1-")

  Extended.neg<-c(T,T,T,F,F,F,F,T,T,T,T,T,T,T,T)

  adductsFile.pos<-data.frame(
    Ion.name=Ion.name.pos[!(Extended.pos & primary)],
    Ion.mass=Ion.mass.pos[!(Extended.pos & primary)],
    Charge=Charge.pos[!(Extended.pos & primary)])
  adductsFile.neg<-data.frame(
    Ion.name=Ion.name.neg[!(Extended.neg & primary)],
    Ion.mass=Ion.mass.neg[!(Extended.neg & primary)],
    Charge=Charge.neg[!(Extended.neg & primary)])

  if(tolower(mode)%in%c("pos","positive","p"))
  {
    # select the adducts (primary or extended)
    tmpAdduct<-adductsFile.pos
    if(!is.na(charge) & is.numeric(charge))
    {
      chargeSign<-NA
      chargeSign<-"+"
      tmpAdduct<-tmpAdduct[tmpAdduct[,"Charge"]==paste(charge,chargeSign,sep=""),]
    }
    if(!is.na(adduct))
    {

      tmpAdduct<-tmpAdduct[tmpAdduct[,"Ion.name"]%in%adduct,]
    }

    if(nrow(tmpAdduct)==0){
      warning("Not found! Reporting all possible adducts")
      tmpAdduct<-adductsFile.pos
    }
    signs<-sapply(stringr::str_extract(as.character(tmpAdduct[,"Ion.mass"]),"\\+|\\-"),function(x){x[[1]]})
    signs[signs=="+"]=-1
    signs[signs=="-"]=+1
    mzCoefficients<-
      stringr::str_extract(sapply(strsplit(as.character(tmpAdduct[,"Ion.mass"]),split = "\\+|\\-"),function(x){x[1]}),
                  "\\d+")

    mzCoefficients[is.na(mzCoefficients)]<-1
    mzCoefficients<-as.numeric(mzCoefficients)
    mzCoefficients[!grepl("/",as.character(tmpAdduct[,"Ion.mass"]),fixed=T)]<-
      1/ mzCoefficients[!grepl("/",as.character(tmpAdduct[,"Ion.mass"]),fixed=T)]

    result<- mz+  (as.numeric(signs)*
                     as.numeric(sapply(strsplit(as.character(tmpAdduct[,"Ion.mass"]),"\\+|\\-"),function(x){x[[2]]})))

    result<- result*mzCoefficients
  }else if(tolower(mode)%in%c("neg","negative","n"))
  {
    tmpAdduct<-adductsFile.neg
    if(!is.na(charge) & is.numeric(charge))
    {
      chargeSign<-NA
      chargeSign<-"-"
      tmpAdduct<-tmpAdduct[tmpAdduct[,"Charge"]==paste(charge,chargeSign,sep=""),]
    }
    if(!is.na(adduct))
    {

      tmpAdduct<-tmpAdduct[tmpAdduct[,"Ion.name"]%in%adduct,]
    }

    if(nrow(tmpAdduct)==0){
      warning("Not found! Reporting all possible adducts")
      tmpAdduct<-adductsFile.neg
    }
    signs<-sapply(str_extract(as.character(tmpAdduct[,"Ion.mass"]),"\\+|\\-"),function(x){x[[1]]})
    signs[signs=="+"]=-1
    signs[signs=="-"]=+1
    mzCoefficients<-
      str_extract(sapply(strsplit(as.character(tmpAdduct[,"Ion.mass"]),split = "\\+|\\-"),function(x){x[1]}),
                  "\\d+")

    mzCoefficients[is.na(mzCoefficients)]<-1
    mzCoefficients<-as.numeric(mzCoefficients)
    mzCoefficients[!grepl("/",as.character(tmpAdduct[,"Ion.mass"]),fixed=T)]<-
      1/ mzCoefficients[!grepl("/",as.character(tmpAdduct[,"Ion.mass"]),fixed=T)]

    result<- mz+  (as.numeric(signs)*
                     as.numeric(sapply(strsplit(as.character(tmpAdduct[,"Ion.mass"]),"\\+|\\-"),function(x){x[[2]]})))
    result<- result*mzCoefficients

  }else
  {
    stop("Incorrect mode! Mode has to be either positive or negative!")
  }


  return(data.frame(correctedMS=result,adductName=tmpAdduct[,"Ion.name"]))

}

#' Calculate back the mz difference given ppm
#'
#' Given mz, charge, ionization, adduct information this fucntion will estimate neutral mass for addcuts
#'
#' @param run A mz value
#' @param ppm ppm difference
#' @return mz difference
#' @keywords internal
ppmCal<-function(run,ppm)
{
  return((run*ppm)/1000000)
}

#' map MSMS data onto camera object
#'
#' This function accepts a CAMERA object, MSMS data and map MSMS data onto the features
#'
#' @param cameraObject A CAMERA object including addcuts and isotopes
#' @param MSMSdata a list of MSnbase objects
#' @param PlusTime The higher boundary of the retention time window
#' @param MinusTime The lower boundary of the retention time window
#' @param ppm ppm mass devation
#' @param listOfMS2Mapped An empty or full list of previously mapped MS2s (the new mapped data will be added to this)
#' @param listOfMS2Mapped An empty or full list of previously unmapped MS2s (the new mapped data will be added to this)
#' @return A list of mapped and/or unmapped MS2s
#' @keywords internal
IntervalMerge<-function(cameraObject,MSMSdata, PlusTime,MinusTime,ppm,listOfMS2Mapped=list(),listOfUnMapped=list()){

  listofPrecursorsmz<-c()
  for(i in seq(1,length(MSMSdata)))
  {
    listofPrecursorsmz<-c(listofPrecursorsmz,MSMSdata[[i]]@precursorMz)
  }

  listofPrecursorsrt<-c()
  for(i in seq(1,length(MSMSdata)))
  {
    listofPrecursorsrt<-c(listofPrecursorsrt,MSMSdata[[i]]@rt)
  }

  CameramzColumnIndex<-which(colnames(cameraObject@groupInfo)=="mz")

  MassRun1<-intervals::Intervals_full(cbind(listofPrecursorsmz,listofPrecursorsmz))

  MassRun2<-intervals::Intervals_full(cbind(cameraObject@groupInfo[,CameramzColumnIndex]-
                                              ppmCal(cameraObject@groupInfo[,CameramzColumnIndex],ppm),
                                            cameraObject@groupInfo[,CameramzColumnIndex]+
                                              ppmCal(cameraObject@groupInfo[,CameramzColumnIndex],ppm)))

  Mass_iii <- intervals::interval_overlap(MassRun1,MassRun2)

  CamerartLowColumnIndex<-which(colnames(cameraObject@groupInfo)=="rtmin")
  CamerartHighColumnIndex<-which(colnames(cameraObject@groupInfo)=="rtmax")

  TimeRun1<-intervals::Intervals_full(cbind(listofPrecursorsrt,listofPrecursorsrt))

  TimeRun2<-intervals::Intervals_full(cbind(cameraObject@groupInfo[,CamerartLowColumnIndex]-MinusTime,
                                            cameraObject@groupInfo[,CamerartHighColumnIndex]+PlusTime))
  Time_ii <- intervals::interval_overlap(TimeRun1,TimeRun2)

  imatch = mapply(intersect,Time_ii,Mass_iii)

  for (i in 1:length(imatch)) {
    for(j in imatch[[i]])
    {
      MSMStmpObject<-MSMSdata[[i]]
      attributes(MSMStmpObject)$fileName<-attributes(MSMSdata)$fileName
      listOfMS2Mapped[[as.character(j)]]<-
        c(listOfMS2Mapped[[as.character(j)]],MSMStmpObject)
    }
  }
  for (i in 1:length(imatch)) {

    if(length(imatch[[i]])==0)
    {
      MSMStmpObject<-MSMSdata[[i]]
      attributes(MSMStmpObject)$fileName<-attributes(MSMSdata)$fileName
      listOfUnMapped<-c(listOfUnMapped,MSMStmpObject)
    }
  }

  return(list(mapped=listOfMS2Mapped,unmapped=listOfUnMapped))
}


#' Map identification results onto camera object
#'
#' This function accepts a CAMERA object, identification results (a data.frame) and map MSMS data onto the features
#'
#' @param metFragSearchResult A data.frame that much contain the identification results and has to have three addional columns: parentMZ, parentRT and fileName
#' @param cameraObject A CAMERA object including addcuts and isotopes
#' @param ppm ppm mass devation
#' @param PlusTime The higher boundary of the retention time window
#' @param MinusTime The lower boundary of the retention time window
#' @param method The method for mapping, This can be one of regular, fast, or parallel
#' @param ncore Number of cores if parallel is selected for mapping method
#' @return A list of Identification where name of each element correspond to an index of a row in cameraObject
#' @keywords internal
#' @import intervals
#' @import CAMERA
#' @import MSnbase
#' @import parallel
map_id_to_camera<-function(metFragSearchResult=NA,cameraObject=NA,ppm=5,MinusTime=5,PlusTime=5,method="regular",ncore)
{
  #metFragSearchResult<-bb
  IDResults<-metFragSearchResult
  #cameraObject<-an
  listofPrecursorsmz<-c()
  listofPrecursorsmz<-IDResults[,"parentMZ"]
  listofPrecursorsrt<-IDResults[,"parentRT"]
  CamerartLowColumnIndex<-which(colnames(cameraObject@groupInfo)=="rtmin")
  CamerartHighColumnIndex<-which(colnames(cameraObject@groupInfo)=="rtmax")
  CameramzColumnIndex<-which(colnames(cameraObject@groupInfo)=="mz")
  imatch=NA
  if(method=="regular")
  {
    MassRun1<-intervals::Intervals_full(cbind(listofPrecursorsmz,listofPrecursorsmz))

    MassRun2<-intervals::Intervals_full(cbind(cameraObject@groupInfo[,CameramzColumnIndex]-
                                     ppmCal(cameraObject@groupInfo[,CameramzColumnIndex],ppm),
                                   cameraObject@groupInfo[,CameramzColumnIndex]+
                                     ppmCal(cameraObject@groupInfo[,CameramzColumnIndex],ppm)))

    Mass_iii <- intervals::interval_overlap(MassRun1,MassRun2)



    TimeRun1<-intervals::Intervals_full(cbind(listofPrecursorsrt,listofPrecursorsrt))

    TimeRun2<-intervals::Intervals_full(cbind(cameraObject@groupInfo[,CamerartLowColumnIndex]-MinusTime,
                                   cameraObject@groupInfo[,CamerartHighColumnIndex]+PlusTime))
    Time_ii <- intervals::interval_overlap(TimeRun1,TimeRun2)

    imatch = mapply(intersect,Time_ii,Mass_iii)
  }else if(method=="fast")
  {
    featureMzs<-cbind(cameraObject@groupInfo[,CameramzColumnIndex]-
                        ppmCal(cameraObject@groupInfo[,CameramzColumnIndex],ppm),
                      cameraObject@groupInfo[,CameramzColumnIndex]+
                        ppmCal(cameraObject@groupInfo[,CameramzColumnIndex],ppm))

    featureRTs<-cbind(cameraObject@groupInfo[,CamerartLowColumnIndex]-MinusTime,
                      cameraObject@groupInfo[,CamerartHighColumnIndex]+PlusTime)

    imatch<-list()
    for(i in 1:length(listofPrecursorsmz))
    {
      mz<-listofPrecursorsmz[i]
      rt<-listofPrecursorsrt[i]

      imatch[[i]]<-which(featureMzs[,1]<mz & featureMzs[,2]>mz & featureRTs[,1]<rt & featureRTs[,2]>rt)


    }
  }else if (method=="parallel"){
    featureMzs<-cbind(cameraObject@groupInfo[,CameramzColumnIndex]-
                        ppmCal(cameraObject@groupInfo[,CameramzColumnIndex],ppm),
                      cameraObject@groupInfo[,CameramzColumnIndex]+
                        ppmCal(cameraObject@groupInfo[,CameramzColumnIndex],ppm))

    featureRTs<-cbind(cameraObject@groupInfo[,CamerartLowColumnIndex]-MinusTime,
                      cameraObject@groupInfo[,CamerartHighColumnIndex]+PlusTime)
    imatch <-parallel::mclapply(c(1:length(listofPrecursorsmz)),FUN =function(x) { which(featureMzs[,1]<listofPrecursorsmz[x] & featureMzs[,2]>listofPrecursorsmz[x] & featureRTs[,1]<listofPrecursorsrt[x] & featureRTs[,2]>listofPrecursorsrt[x])},mc.cores=ncore)

  }else{stop("Method should be either fast,par and regular")}

  listOfMS2Mapped<-list()
  for (i in 1:length(imatch)) {
    for(j in imatch[[i]])
    {
      if(is.null(listOfMS2Mapped[[as.character(j)]]))
      {
        listOfMS2Mapped[[as.character(j)]]<-data.frame(IDResults[i,],stringsAsFactors = F)
      }else
      {
        listOfMS2Mapped[[as.character(j)]]<-
          rbind(listOfMS2Mapped[[as.character(j)]],data.frame(IDResults[i,],stringsAsFactors = F))
      }

    }
  }
  return(list(mapped=listOfMS2Mapped))
}


#' Merge MS2 fragment ions
#'
#' This function accepts a list of MS2s and cluster the corresponding ions based on user defined boundaries.
#'
#' @param spectra A list of MS2s
#' @param eppm ppm differences betweem MS2 ions
#' @param eabs Absolute mass deviation between MS2 ions
#' @param int.threshold The peaks not reaching this intensity will be removed
#' @return A list of merged MS2s similar format as input
#' @keywords internal
#' @import intervals
#' @import CAMERA
#' @import MSnbase
#' @import xcms
merge.spectra.group <- function(spectra, eppm, eabs, int.threshold) {
  if(is.na(spectra) || length(spectra) == 0) return(NULL)
  max.precursor.int <- -1
  rt.at.max.int <- NA
  start.rt <- attributes(spectra[[1]])$rt
  end.rt <- attributes(spectra[[1]])$rt
  all_unique_files<-paste(sapply(spectra,function(x){attributes(x)$fileName}),collapse = ";")
  mean.precursor.mz <- 0
  max.tic <- -1
  tandemms <- lapply(1:length(spectra), function(index) {
    # print(index)
    x<-spectra[[index]]
    # assign local variables
    cur.rt <- attributes(x)$rt
    cur.precursor.int <- attributes(x)$precursorIntensity
    cur.tic <- attributes(x)$tic
    cur.precursor.mz <- attributes(x)$precursorMz
    # assign global variables
    mean.precursor.mz <<- mean.precursor.mz + cur.precursor.mz
    max.tic <<- if(max.tic < cur.tic) cur.tic else max.tic
    if(max.precursor.int < cur.precursor.int) {
      max.precursor.int <<- cur.precursor.int
      rt.at.max.int <<- cur.rt
    }
    start.rt <<- if(start.rt > cur.rt) cur.rt else start.rt
    end.rt <<- if(end.rt < cur.rt) cur.rt else end.rt
    # start merging ms2 data
    a <- cbind(attributes(x)$mz, attributes(x)$intensity)
    colnames(a) <- c("mz", "intensity")
    return(a)
  })
  mean.precursor.mz <- mean.precursor.mz / length(spectra)

  # Filter intensity
  tandemms <- lapply(tandemms, function(spec) spec[spec[,"intensity"]>=int.threshold,])
  tandemrm <- sapply(1:length(tandemms), function(x) { if ( (dim(tandemms[[x]])[1] == 0) || is.null(dim(tandemms[[x]])[1]) ) x <- TRUE else x <- FALSE } )

  if ( (all(as.logical(tandemrm)) == TRUE) | (is.na(all(as.logical(tandemrm)))) )
    return(NULL)
  else
    tandemms <- tandemms[!tandemrm]

  # Process spectra
  peaks <- do.call(rbind, tandemms)
  g <- xcms:::mzClust_hclust(peaks[,"mz"], eppm=eppm, eabs=eabs)

  mz <- tapply (peaks[,"mz"], as.factor(g), mean)
  intensity <- tapply (peaks[,"intensity"], as.factor(g), max)

  if(length(mz) == 0)return(NULL);

  intensity <- intensity[order(mz)]
  mz <- mz[order(mz)]

  res <- .Call("Spectrum2_constructor_mz_sorted",
               attributes(spectra[[1]])$msLevel, length(mz), rt.at.max.int,
               attributes(spectra[[1]])$acquisitionNum,
               attributes(spectra[[1]])$scanIndex, max.tic, mz,
               intensity, attributes(spectra[[1]])$fromFile, attributes(spectra[[1]])$centroided,
               attributes(spectra[[1]])$smoothed, attributes(spectra[[1]])$polarity,
               attributes(spectra[[1]])$merged, attributes(spectra[[1]])$precScanNum,
               mean.precursor.mz, max.precursor.int, attributes(spectra[[1]])$precursorCharge, attributes(spectra[[1]])$collisionEnergy,
               TRUE, attributes(spectra[[1]])$.__classVersion__,
               PACKAGE = "MSnbase")
  attributes(res)$start.rt <- start.rt
  attributes(res)$end.rt <- end.rt
  attributes(res)$fileName <- all_unique_files

  return(res)
}


# filter.grouped.spectra <- function(grouped.spectra, max.rt.range, max.mz.range, min.rt, max.rt, min.mz, max.mz) {
#   filtered.grouped.spectra <- list()
#   sapply(grouped.spectra, function(grouped.spectrum) {
#     mzs<-sapply(grouped.spectrum, function(x) attributes(x)$precursorMz)
#     rts<-sapply(grouped.spectrum, function(x) attributes(x)$rt)
#     mean.mzs <- mean(mzs)
#     mean.rts <- mean(rts)
#     max.number.peaks <- max(sapply(grouped.spectrum, function(x) attributes(x)$peaksCount))
#     if((max(mzs)-min(mzs)) <= max.mz.range & (max(rts)-min(rts)) <= max.rt.range & mean.rts < max.rt & mean.mzs < max.mz
#        & mean.rts > min.rt & mean.mzs > min.mz & max.number.peaks > 0) {
#       filtered.grouped.spectra[[length(filtered.grouped.spectra) + 1]] <<- grouped.spectrum
#     }
#   })
#   return(filtered.grouped.spectra)
# }

#' Calculate back the mz difference given ppm
#'
#' Given mz, charge, ionization, adduct information this fucntion will estimate neutral mass for addcuts
#'
#' @param peak A mz value
#' @param ppm ppm difference
#' @return mz difference
#' @keywords internal
get.ppm.from.abs <- function(peak, ppm) {
  return ((peak / 1000000.0) * ppm)
}

#' Cluster MS2 parent ions
#'
#' This function accepts a list of MS2s and cluster the corresponding parent ions ions based on user defined boundaries.
#'
#' @param spectra A list of MS2s
#' @param mzppm ppm differences betweem MS2 ions
#' @param mzabs Absolute mass deviation between MS2 ions
#' @param rtabs Retention time tolerance
#' @param int.threshold The peaks not reaching this intensity will be removed
#' @return A nested list in which the first level indicate the group of MS2s
#' @keywords internal
#' @import intervals
#' @import CAMERA
#' @import MSnbase
#' @import xcms
collect.spectra.lists <- function(spectra, mzabs, mzppm, rtabs) {
  a<-t(sapply(1:length(spectra), function(spectrum.id) {
    return(c(attributes(spectra[[spectrum.id]])$rt, attributes(spectra[[spectrum.id]])$precursorMz, spectrum.id))
  }))
  a <- a[order(a[,2], a[,1]),]
  index <- 1
  grouped.spectra.tmp <- list()
  while(index <= dim(a)[1]) {
    spectra.list <- list()
    spectra.list[[1]] <- spectra[[a[index,3]]]
    cur.abs <- get.ppm.from.abs(a[index, 2], mzppm * 0.5) + mzabs
    index <- index + 1
    while(index <= dim(a)[1]) {
      diff.mz <- abs(a[index - 1, 2] - a[index, 2])
      if(diff.mz <= cur.abs) {
        spectra.list[[length(spectra.list) + 1]] <- spectra[[a[index,3]]]
        index <- index + 1
      } else break
    }
    grouped.spectra.tmp[[length(grouped.spectra.tmp) + 1]] <- spectra.list
  }
  grouped.spectra <- list()
  for(spectrum.group.index in 1:length(grouped.spectra.tmp)) {
    spectrum.group <- grouped.spectra.tmp[[spectrum.group.index]]
    a<-t(sapply(1:length(spectrum.group), function(spectrum.index) {
      c(attributes(spectrum.group[[spectrum.index]])$rt, attributes(spectrum.group[[spectrum.index]])$precursorMz, spectrum.index)
    }))
    a <- matrix(a[order(a[,1], a[,2]),], ncol=3)
    index <- 1
    while(index <= dim(a)[1]) {
      spectra.list <- list()
      spectra.list[[1]] <- spectrum.group[[a[index,3]]]
      index <- index + 1
      while(index <= dim(a)[1]) {
        diff.rt <- abs(a[index - 1, 1] - a[index, 1])
        if(diff.rt <= rtabs) {
          spectra.list[[length(spectra.list) + 1]] <- spectrum.group[[a[index,3]]]
          index <- index + 1
        } else break
      }
      grouped.spectra[[length(grouped.spectra) + 1]] <- spectra.list
    }
  }
  return(grouped.spectra)
}

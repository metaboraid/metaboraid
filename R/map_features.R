#' Mapping MS2 to CAMERA
#'
#' This function allows you map MS2 data onto the CAMERA quantification data.
#' We use this function in downstream to create HyperMS2 objects and map data to addcuts etc.
#'
#' @param inputMS2s A character vector containing the absolute paths to the MS2 mzML files.
#' @param camera_results A CAMERA object e.g. from running CAMERA on the data
#' @param ppm ppm tolerance between the parent ion and the center of the features. Default 10
#' @param rt retention time tolerance between the parent ion and the center of the features in seconds. Default 10
#' @param verbose Show information about different stages of the processes. Default FALSE
#' @export
#' @examples
#'
#' library(CAMERA)
#' library(metaboraid)
#' # Read MS1 and MS2 files
#' ms1_files<-system.file("ms1data",c("X1_Rep1.mzML","X2_Rep1.mzML"),package = "metaboraid")
#' ms2_files<-system.file("ms2data",c("sample1.mzML","sample2.mzML"),package = "metaboraid")
#'
#' # mass trace detection
#' xs <- xcmsSet(ms1_files,method="centWave",ppm=30,peakwidth=c(5,10))
#'
#' # mass trace matching
#' xsg <- group(xs)
#'
#' # convert to CAMERA
#' xsa <- xsAnnotate(xsg)
#'
#' # Group mass traces
#' anF <- groupFWHM(xsa, perfwhm = 0.6)
#'
#' # Detect isotopes
#' anI <- findIsotopes(anF, mzabs = 0.01)
#'
#' # Group using correlation
#' anIC <- groupCorr(anI, cor_eic_th = 0.75)
#'
#' # Find adducts
#' anFA <- findAdducts(anIC, polarity="positive")
#'
#' # map features and MS2s
#' mapped_features<-map_features(inputMS2s = ms2_files,input_camera = anFA,ppm = 10,rt = 10)
#'

#' @return
#' A list with two named elements (mapped, unmapped). The mapped contain a nested list where each element has a name and contain a list of MS2 spectra for a CAMERA feature with index==name of the list.
#' The unmapped element contains the crude list of all the unmapped MS2s
#'
#' @import CAMERA
#' @import xcms
#' @import MSnbase
#'

map_features<-function(inputMS2s=NA,input_camera=NA,ppm=10,rt=10,
                       verbose=T)
{
  if(any(is.na(inputMS2s)) | any(is.null(inputMS2s)))
  {
    stop("No MS2 files have been provided!")
  }

  if(any(is.na(input_camera)) | any(is.null(input_camera)))
  {
    stop("No CAMERA input have been provided!")
  }

  if(any(!file.exists(inputMS2s)))
  {
    writeLines(sprintf("\"%s\" was not found!",(inputMS2s[!file.exists(inputMS2s)])))
  }

  if(class(input_camera)[[1]]!="xsAnnotate")
  {
    stop("input_camera is not a CAMERA object!")
  }

  tmp<-CAMERA::getPeaklist(input_camera)
  if(any(!c("isotopes","adduct","pcgroup")%in%colnames(tmp)))
  {
    stop("CAMERA object does not contain pcgroup,isotopes or adduct! Remember to run grouping, isotope and adduct detection using CAMERA before using this function!")

  }
  rm(tmp)

  if(any(is.na(rt)) | any(is.null(rt)))
  {
    stop("No rt input have been provided!")
  }

  if(any(is.na(ppm)) | any(is.null(ppm)))
  {
    stop("No ppm input have been provided!")
  }

  if(!is.numeric(rt))
  {
    stop("rt must be numberic")
  }

  if(!is.numeric(ppm))
  {
    stop("ppm must be numberic")
  }

  if(ppm<0)
  {
    stop("ppm must be numberic and higher than 0")
  }

  if(rt<0)
  {
    stop("rt must be numberic and higher than 0")
  }


  mappingResult<-list(mapped=list(),unmapped=list())

  for(MS2 in inputMS2s)
  {
    if(verbose)cat("\nReading MS2 file: ",MS2,"\n")
    MSMSdata<-MSnbase::readMSData(MS2, msLevel = 2, verbose = verbose)
    attributes(MSMSdata)$fileName<-basename(MS2)
    if(verbose)cat("\nMapping to MS2 to CAMERA: ",MS2,"\n")
    mappingResult<-IntervalMerge(cameraObject = input_camera,MSMSdata = MSMSdata,PlusTime = rt,
                                 MinusTime = rt,ppm = ppm,listOfMS2Mapped = mappingResult$mapped,
                                 listOfUnMapped = mappingResult$unmapped)
  }
  if(verbose)cat("\nNumber of mapped MS2 ions: ",length(mappingResult$mapped),
                 " Number of unmapped MS ions: ",length(mappingResult$unmapped),"\n")
  return(mappingResult)
}

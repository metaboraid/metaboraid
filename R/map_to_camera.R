#' Map the identification results to quantification
#'
#' This function allows you to map the identification results to quantification data from CAMERA.
#'
#'map_to_camera<-function(input_camera=NA,
#' @param input_camera A CAMERA object (results of preprocessing using CAMERA)
#' @param identification_result dataframe for example from running one of the search engines
#' @param ppm Precursor ppm mass deviation for matching features and identification parent ions. Default 10
#' @param rt Retention time tolerance for matching features and identification parent ions. Default 10
#' @param higher_the_better If the score of the search engine is higher the better score. Default TRUE
#' @param score_column Which column of identification_result shows the score.
#' @param impute Whether to impute the unknown mass traces using PC groups
#' @param mapping_method Mapping method to use. It has to be one regular, fast, or parallel.
#' @param ncores Number of cores if mapping_method is parallel
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
#' # Map adducts
#' mapped_adducts<-map_adducts(inputMS2List=mapped_features,input_camera=anFA,
#'                             precursorppm=10,
#'                             fragmentppm=20,fragmentabs=0.01,minPrecursorMass=NA,maxPrecursorMass=NA,
#'                             minPeaks=10,maxSpectra=10,mode="pos",adductRules="primary",
#'                             outputDir="general_parameters_4",searchMultipleChargeAdducts=T,
#'                             includeMapped=T,includeUnmapped=F,verbose=T)
#'
#' # Run the search
#' ID_results<-run_sirius("parameter_files.zip",database = "all",ncores = 2,progress_bar = T,verbose = T,results_folder = "results",chech_file_interval = 2,timeout = 600,conda = "auto")
#'
#' # map ID to features
#'
#' final_results<-metaboraid:::map_to_camera(input_camera = anFA,identification_result = ID_results,ppm = 10,rt = 10,higher_the_better = T,score_column = "SiriusScore",impute = F,mapping_method = "fast")
#'
#' @return
#' A list with two elements: peakMatrix, variableData. peakMatrix contains the peak information from CAMERA and variableData includes the identification information. The first columns of the both matrices are the key to match these two dataframes.
#'
#' @details
#' This is normally the last function applied on the data. The purpose of this function is to provide ready to analyze data by simply finding the identified features.
#' The option impute uses CAMERA PC groups to impute identification results within this group. For example if one of the members has been identified, the rest will be imputed based on this member.
#' We always take the top scoring hit for imputing and for mapping
#'
#' @import reticulate
#' @import zip
#' @import parallel
#' @import future
#' @import progressr
#' @import future.apply
#'

map_to_camera<-function(input_camera=NA,
                        identification_result=NA,ppm=10,rt=10,higher_the_better=TRUE,score_column=NA,impute=FALSE,mapping_method="fast",ncores=1,verbose=FALSE){

  if(verbose)cat("Checking inputs ...","\n")
  if(is.na(score_column)|is.null(score_column))
  {
    stop("score_column has not been provided")
  }

  if(!is.matrix(identification_result) & !is.data.frame(identification_result))
  {
    stop("No identification_result has been provided!")
  }else if(nrow(identification_result)<1){
    stop("Empty IDs in identification_result!")
  }else if(!all(c("parentMZ","parentRT","fileName",score_column)%in%colnames(identification_result))){
    stop(sprintf("Identification results must contain at least these columns: %s!",paste(c("parentMZ","parentRT","fileName",score_column),collapse = ", ")))
  }

  if(any(is.na(input_camera)) | any(is.null(input_camera)))
  {
    stop("No input_camera input have been provided!")
  }


  if(class(input_camera)[[1]]!="xsAnnotate")
  {
    stop("input_camera is not a CAMERA object!")
  }

  cameraPeakList<-CAMERA::getPeaklist(input_camera)
  #if(any(!c("mz","mzmin","mzmax","rt","rtmin","rtmax","npeaks","isotopes","adduct","pcgroup")%in%colnames(cameraPeakList)))
  #{
  #  stop("CAMERA object must contain at least mz, mzmin, mzmax, rt, rtmin, rtmax, npeaks, isotopes, adduct, and pcgroup! Remember to run grouping, isotope and adduct detection using CAMERA before using this function!")

  #}


  if(is.na(mapping_method)|is.null(mapping_method))
  {
    stop("mapping_method must be one of regular, fast or parallel!")

  }else if((!any(mapping_method%in%c("regular","fast","parallel"))) | length(mapping_method)>1)
  {
    stop("mapping_method must be one of regular, fast or parallel!")
  }

  if(is.na(ncores)|is.null(ncores))
  {
    stop("ncores must be a number greater than 0!")
  }else if(ncores<1){
    stop("ncores must be equal or greater than 0!")
  }


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

  if(is.null(higher_the_better)|is.na(higher_the_better)){
    stop("higher_the_better must be either TRUE or FALSE!")
  }else if(!is.logical(higher_the_better)){
    stop("higher_the_better must be either TRUE or FALSE!")
  }

  if(verbose)cat("Extracting peak list ...","\n")
  cameraInformation<-cameraPeakList[,c("mz","mzmin","mzmax","rt","rtmin","rtmax","npeaks","isotopes","adduct","pcgroup")]
  colnames(cameraInformation)<-paste("xcmsCamera_",colnames(cameraInformation),sep = "")


  if(verbose)cat("Mapping ...","\n")
  mappedToCamera<-map_id_to_camera(metFragSearchResult = identification_result,
                                  cameraObject = input_camera,MinusTime = rt,PlusTime = rt,ppm = ppm,method=mapping_method,ncore = ncores)


  VariableData<-data.frame(matrix("Unknown",nrow = nrow(cameraPeakList),
                                  ncol = (ncol(identification_result)+1+ncol(cameraInformation))),stringsAsFactors = F)

  colnames(VariableData)<-c("variableMetadata",colnames(identification_result),colnames(cameraInformation))

  VariableData[,"variableMetadata"]<-paste("variable_",1:nrow(cameraPeakList),sep="")
  VariableData[,colnames(cameraInformation)]<-cameraInformation

  for(rnName in rownames(cameraPeakList))
  {
    if(rnName %in% names(mappedToCamera$mapped))
    {
      tmpId<-mappedToCamera$mapped[[rnName]]
      if(higher_the_better)
      {
        tmpId<-tmpId[which.max(tmpId[,score_column]),]
      }else{
        tmpId<-tmpId[which.min(tmpId[,score_column]),]
      }

      VariableData[VariableData[,"variableMetadata"]==paste("variable_",rnName,sep=""),
                   c(2:(ncol(identification_result)+1))]<-tmpId

    }
  }


  VariableData$imputed<-"No"

  if(impute)
  {
    if(verbose)cat("imputting ...","\n")
    toBeImputed<-which(VariableData[,2]=="Unknown")
    pcgroups<-VariableData[toBeImputed,"xcmsCamera_pcgroup"]

    for(pcgr in unique(pcgroups))
    {
      selectedFeatures<-
        VariableData[,"variableMetadata"]%in%(VariableData[VariableData[,"xcmsCamera_pcgroup"]==pcgr,"variableMetadata"]) &
        VariableData[,"parentMZ"]!="Unknown"

      if(any(selectedFeatures))
      {
        tmpIDs<-VariableData[selectedFeatures,]
        tmpId<-NA
        if(higher_the_better)
        {
          tmpId<-tmpIDs[which.max(tmpIDs[,score_column]),]
        }else
        {
          tmpId<-tmpIDs[which.min(tmpIDs[,score_column]),]
        }


        imputedVariables<- (VariableData[VariableData[,"xcmsCamera_pcgroup"]==pcgr,"variableMetadata"])

        imputedVariables<- VariableData[,"variableMetadata"]%in%imputedVariables & VariableData[,2]=="Unknown"

        VariableData[imputedVariables,c(2:(ncol(identification_result)+1))]<-tmpId[,c(2:ncol(tmpId))]
        VariableData[imputedVariables,"imputed"]<-"yes"

      }
    }

  }
  if(verbose)cat("Formatting ...","\n")
  peakMatrix<-cbind.data.frame(dataMatrix=VariableData[,"variableMetadata"],cameraPeakList,stringsAsFactors = F)
  VariableData<-sapply(VariableData, gsub, pattern="\'|#", replacement="")
  return(list(peakMatrix=peakMatrix,variableData=VariableData))
}


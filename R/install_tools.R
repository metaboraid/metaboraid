#' Install tools needed by MetaboRAID
#'
#' This function allows you to install the tools needed for MetaboRAID.
#' We use conda to install the tools so make sure you have the latest version conda installed (\url{https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html}).
#'
#' At this stage, CSI-FINGERID, Metfrag and Chemdistiller can be installed on Linux, macOS and Windows however, you should install CFM-ID yourself (\url{https://cfmid.wishartlab.com/})
#'
#' We install everything in metaboraid_package environment. However, CFM-ID will be installed in a different environment.
#'
#' @param method Which method to use to install the packages (only conda is supported)
#' @param conda path to conda binary. Default: auto
#' @param version Version of the package to install
#' @param cfm_id Set to TRUE if you want to run CFM-ID. Please remember that this is not supported on Windows!
#' @export
#' @examples
#' install_tools()
#'
#' @return
#' none
#' @import reticulate
install_tools <- function(method = c("conda"),
                          conda = "auto",
                          version = "default", cfm_id=FALSE,
                          ...) {

  # verify method
  method <- match.arg(method)

  # resolve version
  if (identical(version, "default"))
    version <- ""
  else
    version <- paste0("==", version)

  # some special handling for windows
  if (is_windows()) {

    # conda is the only supported method on windows
    method <- "conda"

    # confirm we actually have conda
    have_conda <- !is.null(tryCatch(conda_binary(conda), error = function(e) NULL))
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
  }

  reticulate::conda_create(
    envname = "metaboraid_package",
    #python_version = "3.8" #worked for Ubuntu in case there are conda package version conflicts
    packages = "python",
    forge = TRUE,
    channel = character(),
    conda = "auto"
  )


  reticulate::conda_install(
    packages = "metaboraid",
    forge = TRUE,
    channel = c("payamemami","bioconda","conda-forge","anaconda","defaults"),envname = "metaboraid_package"
  )

  if(!is.na(cfm_id) & !is.null(cfm_id))
  {

  }else{
    if(cfm_id==TRUE){
      if(is_windows()){
        cat("We cannot install cfm-id using conda on Windows. Please install it manually. Visit https://cfmid.wishartlab.com/")

      }else{

        reticulate::conda_install(
          packages = "cfm",
          forge = TRUE,
          channel = c("bioconda","conda-forge","anaconda","defaults"),envname = "metaboraid_package_cfm"
        )
      }

    }
  }

  check_if_software_exist("metfrag")
  check_if_software_exist("sirius")

}


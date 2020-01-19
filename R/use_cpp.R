#' Function for activating cpp file for samuuika with TMB
#'
#' @import TMB
#' @importFrom TMB compile
#' @importFrom TMB dynlib
#' @param TmbFile Cpp file name
#' @param CppDir directory having the cpp file
#' @param RunDir directory for running the code (default: working directory)
#' @param overwrite whether the cpp file is overwritten (TRUE) or not (FALSE; default)
#'
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#' use_samuika_tmb()
#' }
#'
#' @export

use_saumika_tmb <- function(TmbFile = "samuika",
                         CppDir = system.file("inst/executable",package="messir"),
                         RunDir = getwd(),
                         overwrite = FALSE) {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Please install TMB package!")
  }
  file.copy( from=paste0(CppDir,"/",TmbFile,".cpp"), to=paste0(RunDir,"/",TmbFile,".cpp"), overwrite=overwrite)
  TMB::compile( paste0(TmbFile,".cpp") )
  dyn.load(TMB::dynlib(TmbFile))
}

#' Function for activating cpp file for samuuika with TMB
#'
#' @import TMB
#' @importFrom TMB compile
#' @importFrom TMB dynlib
#' @param TmbFile Cpp file name
#' @param CppDir directory having the cpp file
#' @param RunDir directory for running the code (default: working directory)
#' @param overwrite whether the cpp file is overwritten (TRUE) or not (FALSE; default)
#'
#' @encoding UTF-8
#'
#' @examples
#' \dontrun{
#' use_est_MSY_tmb()
#' }
#'
#' @export

use_est_MSY_tmb <- function(TmbFile = "est_MSY",
                            CppDir = system.file("inst/executable",package="messir"),
                            RunDir = getwd(),
                            overwrite = FALSE) {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Please install TMB package!")
  }
  file.copy( from=paste0(CppDir,"/",TmbFile,".cpp"), to=paste0(RunDir,"/",TmbFile,".cpp"), overwrite=overwrite)
  TMB::compile( paste0(TmbFile,".cpp") )
  dyn.load(TMB::dynlib(TmbFile))
}

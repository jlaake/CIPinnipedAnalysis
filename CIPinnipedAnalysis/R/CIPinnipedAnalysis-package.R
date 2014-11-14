#' CIPinnipedAnalysis
#' 
#' Channel Islands Pinniped Analysis functions
#' : pup production, early pup mortality, SST and other enviromental measurement computations
#' 
#' \tabular{ll}{ Package: \tab CIPinnipedAnalysis\cr Type: \tab Package\cr
#' Version: \tab 1.0\cr Date: \tab 2009-11-30\cr License: \tab GPL-2\cr
#' LazyLoad: \tab yes\cr }
#' 
#' @name CIPinnipedAnalysis-package
#' @aliases CIPinnipedAnalysis-package CIPinnipedAnalysis
#' @docType package
#' @author Jeff Laake
#' 
#' Maintainer: <jeff.laake@@noaa.gov>
#' @keywords package
#' @examples
#' # Note this uses default argument value of NULL for dir which
#' # means that the locations for the databases is specified by databases.txt in CalcurData package
#' check.dead.pup.areas()
#' prod=do.pup.production()
#' do.early.pup.mortality()
#' 

NULL

#' POPAN model results and correction factors from dead pup tagging data
#' 
#' POPAN results and correction factor data (estimates) that are used to correct dead
#' pup counts for years without tagging data.
#' 
#' @name smi1994.popan.results
#' @aliases smi1994.popan.results smi1995.popan.results smi1998.popan.results smi2002.popan.results sni2006.popan.results
#' @docType data
#' @format A list with elements: popan.results - a marklist of fitted POPAN models, cfdata - list of correction factor data used to correct
#' dead pup counts for years without tagging data, and cfbyocc - dataframe of correction factors for each occasion. 
#' @keywords datasets

NULL






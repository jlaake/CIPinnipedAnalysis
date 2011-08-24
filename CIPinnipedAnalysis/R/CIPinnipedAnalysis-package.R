

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
#' 
#' # fdir should be set to location of  BrandMaster.mdb
#' # each defaults below to current working directory
#' # Use "J:/" if it should be Calcurr/Databases
#' # fdir="J:/"
#' fdir="" 
#' getwd()
#' check.dead.pup.areas(fdir=fdir)
#' do.pup.production(fdir=fdir)
#' do.early.pup.mortality(fdir=fdir)
#' example(create.SST.anomalies)
#' 
NULL




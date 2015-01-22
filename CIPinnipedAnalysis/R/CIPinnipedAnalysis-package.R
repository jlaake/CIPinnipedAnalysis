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
#' @keywords package
#' @import RMark CalcurData plotrix nlme
#' @examples
#' # Note this uses default argument value of NULL for dir which
#' # means that the locations for the databases is specified by databases.txt in CalcurData package
#' check.dead.pup.areas()
#' source(file.path(system.file(package="CIPinnipedAnalysis"),"PupProduction.r"))
#' source(file.path(system.file(package="CIPinnipedAnalysis"),"EarlyPupMortality.r"))
#' 

NULL

#' POPAN model results and correction factors from dead pup tagging data
#' 
#' POPAN results and correction factor data (estimates) that are used to correct dead
#' pup counts for years without tagging data.
#' 
#' @name smi1994.popan.results
#' @aliases smi1994.popan.results smi1995.popan.results smi1998.popan.results smi2002.popan.results smi1998a.popan.results smi2002a.popan.results sni2006.popan.results
#' @docType data
#' @format A list with elements: popan.results - a marklist of fitted POPAN models, cfdata - list of correction factor data used to correct
#' dead pup counts for years without tagging data, and cfbyocc - dataframe of correction factors for each occasion.  The smi1998a.popan.results and smi2002a.popan.results
#' use areas for predictions instead of covariates to enable corrections for years <1997.
#' @keywords datasets

NULL


#' Calcofi Environmental Data Measurements
#' 
#' A sequence of oceanographic measurements at seven calcofi stations in July from 1984-2013
#' Data were provided by Isaac Schroeder.
#' 
#' @name calcofi
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' Month - either July or October
#' Year - year of data measurement
#' Station - station value
#' dynamic_height_0_500m - dynamic heigth measurement - units unknown; function of temp and salinity
#' stratification 
#' pycnocline_depth  
#' R_POTEMP_25m - temperature at 25 meters
#' R_POTEMP_75m - temperature at 75 meters
#' R_SIGMA_25m - density gradient at 25m
#' R_SIGMA_75m - density gradient at 75m
#' R_O2_25m - oxygen at 25m
#' R_O2_75m - oxygen at 75m
#' R_NO3_25m - NO3 at 25m
#' R_NO3_75m - NO3 at 75m
#' R_SIGMA_75m - density gradient at 75m}
#' @keywords datasets

NULL



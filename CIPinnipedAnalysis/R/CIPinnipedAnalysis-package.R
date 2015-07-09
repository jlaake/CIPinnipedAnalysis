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

#' Port San Luis Sea Level Height Data
#' 
#' A sequence of monthly average sea level height at Port San Luis
#' http://tidesandcurrents.noaa.gov/sltrends/sltrends_station.shtml?stnid=9412110 
#' 
#' @name SeaLevelHeight
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' Year - year of data measurement
#' Month - numeric month of measurement
#' SeaLevelHeight - monthly mean sea level height in meters with seasonal fluctuations removed
#' } 
#' @keywords datasets

NULL



#' Scat record data
#' 
#' A record for each scat sampled 
#' 
#' @name scats
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' SCATNUM - scat identifier
#' SPCOD   - code for pinniped species
#' DATE    - date scat collected
#' OTO     - yes of TRUE if any otloliths in scat
#' BKS     - yes or TRUE if any beaks in scat
#' BNE     - yes or TRUE if any bone in scat
#' Year     - year of collection
#' Month    - month of collection
#' }
#' Created with scats.txt file and the collowing R code
#' scats=read.delim("scats.txt",sep="\\t",header=TRUE,colClasses=c("factor","factor","character","NULL","NULL","NULL","NULL","character","character","character","numeric","numeric"))
#' #remove 1994 and 2007 because there is no prey species fo data
#' scats=scats[!scats$Year%in%c(1994,2007),]
#' # exclude any with all FALSE or NO for OTO,BKS,BNE except those in include which have prey data
#' xx=merge(scats,fo)
#' include=xx[(toupper(xx$BNE)%in%c("FALSE","NO")&toupper(xx$OTO)%in%c("FALSE","NO")&toupper(xx$BKS)%in%c("FALSE","NO"))&!is.na(xx$PREYSP),"SCATNUM"]
#' rm(xx)
#' scats=scats[toupper(scats$BNE)%in%c("YES","TRUE")|toupper(scats$BKS)%in%c("YES","TRUE")|toupper(scats$OTO)%in%c("YES","TRUE")|scats$SCATNUM%in%include,]
#' scats=droplevels(scats)
#' save(scats,file="scats.rda")
#' @keywords datasets

NULL

#' Prey species ID data
#' 
#' Prey species identified from Zalophus scat data
#' 
#' @name fo
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' SCATNUM - scat identifier
#' PREYSP  - code for prey species in scat
#' FO_OT.BK - 1 if identified by beak or otolith
#' FO_OTHER - 1 if indentified by other structure
#' Date     - date scat collected
#' Year     - year of collection
#' }
#' Created with FO_DAT.txt file and following R code
#' fo=read.delim("Fo_DAT.txt",sep="\\t",header=TRUE,colClasses=c("factor","factor","character","character","NULL","character","NULL","NULL","numeric","NULL"))
#' #remove 1994 and 2007 because there is no prey species fo data
#' fo=droplevels(fo[!fo$Year%in%c(1994,2007),])
#' save(fo,file="fo.rda")
#' @keywords datasets

NULL


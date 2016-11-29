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
#' @importFrom graphics abline lines plot
#' @importFrom stats na.exclude AIC approx coef cov prcomp predict var
#' @importFrom utils data
#' @examples
#' # Note this uses default argument value of NULL for dir which
#' # means that the locations for the databases is specified by databases.txt in CalcurData package
#' \donttest{
#' check.dead.pup.areas()
#' source(file.path(system.file(package="CIPinnipedAnalysis"),"PupProduction.r"))
#' source(file.path(system.file(package="CIPinnipedAnalysis"),"EarlyPupMortality.r"))
#' }

NULL

#' Zc production and early mortality
#' 
#' @name Zc
#' @examples
#' \donttest{
#'  check.dead.pup.areas()
#'  source(file.path(system.file(package="CIPinnipedAnalysis"),"ZcPupProduction.r"))
#'  source(file.path(system.file(package="CIPinnipedAnalysis"),"ZcEarlyPupMortality.r"))
#' }

NULL

#' Cu production and early mortailty
#' 
#' @name Cu
#' @examples
#' \donttest{
#'  check.dead.pup.areas()
#'  source(file.path(system.file(package="CIPinnipedAnalysis"),"CuPupProduction.r"))
#'  source(file.path(system.file(package="CIPinnipedAnalysis"),"CuEarlyPupMortality.r"))
#' }

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


#' Anchovy Biomass Estimates
#' 
#' A record for 1972-2011 from
#' MacCall, A. D., W. J. Sydeman, P. C. Davison, and J. A. Thayer. 2016. Recent collapse of northern anchovy biomass off California. Fish. Res. 175:87-94.
#' provided as file from Bill Sydeman
#' 
#' @name anchovy
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' Year           - year 
#' AnchovyBiomass - estimate of stock biomass 1000 metric tons for the year
#' }
#' @keywords datasets

NULL

#' Hake Biomass Estimates
#' 
#' A record for 1977-2015 from
#' Grandin, C.J., A.C. Hicks, A.M. Berger, A.M. Edwards, N. Taylor, I.G. Taylor, and S. Cox. 2016. Status of the Pacific Hake (whiting) stock in U.S. and Canadian waters in 2016. Prepared by the Joint Technical Committee of the U.S. and Canada Pacific Hake/Whiting Agreement, National Marine Fisheries Service and Fisheries and Oceans Canada. 165 p.
#' 
#' @name hake
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' Year           - year 
#' Age0           - biomass of Age 0 hake
#' Age1           - biomass of Age 1 hake
#' Age2           - biomass of Age 2 hake
#' }
#' @keywords datasets

NULL

#' Sardine Biomass Estimates
#' 
#' A record for 1981-2013 of age 0 and age 1 sardines from
#' Hill, K. T., P. R. Crone, E. Dorval, and Macewicz. 2015. ASSESSMENT OF THE PACIFIC SARDINE RESOURCE IN 2015 FOR U.S.A. MANAGEMENT IN 2015-16. NOAA-TM-NMFS-SWFSC-546 .
#' 
#' @name sardine
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' Year           - year 
#' Age            - either age 0 or age 1
#' Abundance      - estimate of abundance at age
#' Biomass        - estimate of stock biomass at age for the year
#' }
#' @keywords datasets

NULL

#' Loligo Quarterly Biomass Estimates
#' 
#' From region 2, using M=0.15 Table 2 from
#' Dorval, E., P. R. Crone, and J. D. McDaniel. 2013. Variability of egg escapement, fishing mortality and spawning population in the market squid fishery in the California Current Ecosystem. Mar. Freshw. Res. 64:80-90.
#' A record for 1999-2006  
#' 
#' @name loligo
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' Year           - year 
#' Quarter        - quarter of year
#' SSBmt          - stock biomass estimate in metric tons
#' }
#' @keywords datasets

NULL

#' Rockfish Biomass Estimates
#' 
#' A record for 1980-2010 from Table 7 of
#' Ralston, S., K. M. Sakuma, and J. C. Field. 2013. Interannual variation in pelagic juvenile rockfish (Sebastes spp.) abundance - going with the flow. Fish. Oceanogr. 22:288-308.
#' 
#' @name Rockfish
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' Year           - year 
#' entomelas      - recruitment estimate/1000 of the species
#' goodei         - recruitment estimate/1000 of the species
#' jordani        - recruitment estimate/1000 of the species
#' paucispinis    - recruitment estimate/1000 of the species
#' pinniger       - recruitment estimate/1000 of the species
#' }
#' @keywords datasets

NULL

#' CALCOFI Fisheries Landings
#' 
#' A record for 1975-2014 
#' 
#' @name Landings
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' Year           - year 
#' Sardine        - landings of sardine
#' Anchovy        - landings of anchovy
#' Mackerel       - landings of mackerel
#' Squid          - landings of squid
#' }
#' @keywords datasets

NULL

#' Callorhinus frequency of occurrence data
#' 
#' A record for each species in a scat
#' 
#' @name cufo
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' SCATNUM           - scat number id 
#' PREYSP            - prey species
#' FO_OT.BK          - landings of anchovy
#' FO_OTHER          - landings of mackerel
#' DATE              - landings of squid'
#' Year              - 
#' }
#' @keywords datasets

NULL

#' Zalophus stranding data from Marine Mammal Center
#' 
#' Provided by Denise Greig from her stranding paper
#' Greig, D.J., F.M.D. Gulland, and C. Kreuder. 2005. A decade of live California sea lion
#'(Zalophus californianus) strandings along the central California coast: Causes and trends,
#' 1991-2000. Aquat. Mammals 31: 11-22. 
#' 
#' @name strandingAL
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' Name           - sea lion identifier 
#' Age            - sea lion age in years
#' Sex            - sex of sea lion M (male) or F (Female) 
#' Length         - standard length in cm
#' }
#' @keywords datasets

NULL


#' Daily Upwelling Index values at 33N
#' http://www.pfeg.noaa.gov/products/PFEL/modeled/indices/upwelling/NA/data_download.html
#' 
#' A record for from 1967 onwards
#' 
#' @name uwi33
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' YYYYMMDD    - 4 digit year, 2 digit month and 2 digit day
#' Index       - upwelling index for the day 
#' }
#' @keywords datasets

NULL


#' Daily Upwelling Index values at 36N
#' http://www.pfeg.noaa.gov/products/PFEL/modeled/indices/upwelling/NA/data_download.html
#' 
#' A record for from 1967 onwards
#' 
#' @name uwi36
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' YYYYMMDD    - 4 digit year, 2 digit month and 2 digit day
#' Index       - upwelling index for the day 
#' }
#' @keywords datasets

NULL


#' Daily Upwelling Index values at 39N
#' http://www.pfeg.noaa.gov/products/PFEL/modeled/indices/upwelling/NA/data_download.html
#' 
#' A record for from 1967 onwards
#' 
#' @name uwi39
#' @docType data
#' @format A dataframe with fields 
#' \preformatted{
#' YYYYMMDD    - 4 digit year, 2 digit month and 2 digit day
#' Index       - upwelling index for the day 
#' }
#' @keywords datasets

NULL


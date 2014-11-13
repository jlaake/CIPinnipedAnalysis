#' Extracts weights for Cu pups
#' : from various tables in the cutagnew.mdb ACCESS
#' database and returns a dataframe with ancillary fields added.
#' 
#' 
#' @param fdir directory for cutagnew.mdb
#' @param ENYears values of years that will be flagged as ENSO years
#' @return dataframe of weights for Cu
#' @export
#' @author Jeff Laake
#' @examples 
#' source(file.path(system.file(package="CIPinnipedAnalysis"),"CuWeightAnalysis.r"))
#' 
get.cu.weights=function(fdir=NULL,ENYears=c(1976,1983,1984,1986,1987,1992,1997,1998,2002,2009))
{
#
# read in weights from cutags table of the Access database
#
cuweights.acv=getCalcurData("Cu","Cutags",dir=fdir)
#
# Set maxyear to largest year; optionally you can set a lower bound
#
maxyear=max(cuweights.acv$cohort)
#
# Exclude Castle Rock and those after maxyear
#
cuweights.acv=cuweights.acv[cuweights.acv$area%in%c("ACV","EAC")&cuweights.acv$cohort<=maxyear,]
#
#  Exclude those with missing weight and unknown sex
#
cuweights.acv=cuweights.acv[!cuweights.acv$sex=="U" & !is.na(cuweights.acv$weight),]
cuweights.acv$sex=factor(cuweights.acv$sex)
#
# Exclude weights before 1 Sept
#
cuweights.acv=cuweights.acv[cuweights.acv$sitedate>as.POSIXct(paste(as.character(cuweights.acv$cohort),"-09-01",sep="")),]
#
#  read in the discard weights from the UnmarkedPupWeights table and exclude any > maxyear
#
cuweights.unmark=getCalcurData("Cu","UnmarkedPupWeights",dir=fdir)
cuweights.unmark=cuweights.unmark[!cuweights.unmark$sex=="U" & !is.na(cuweights.unmark$weight),]
cuweights.unmark=cuweights.unmark[cuweights.unmark$cohort<=maxyear,]
# Read duplicate tag table
cuweights.dup=getCalcurData("Cu","DuplicateTags0708",dir=fdir)
cuweights.dup=cuweights.dup[!cuweights.dup$sex=="U" & !is.na(cuweights.dup$weight),]
#
#   Add a days field which is the number of days from 1 Oct of the year
#
cuweights.acv$days=floor(as.numeric((cuweights.acv$sitedate-as.POSIXct(paste(as.character(cuweights.acv$cohort),"-10-01",sep=""))))/(24*3600))
cuweights.unmark$days=floor(as.numeric((cuweights.unmark$sitedate-as.POSIXct(paste(as.character(cuweights.unmark$cohort),"-10-01",sep=""))))/(24*3600))
cuweights.dup$days=floor(as.numeric((cuweights.dup$sitedate-as.POSIXct(paste(as.character(cuweights.dup$cohort),"-10-01",sep="")))))
#
#  Combine data tables
#
cuweights.dup=subset(cuweights.dup,select=c("cohort","sex","weight","length","girth","pelage","days"))
cuweights.unmark=subset(cuweights.unmark,select=c("cohort","sex","weight","length","girth","pelage","days"))
cuweights=subset(cuweights.acv,select=c("cohort","sex","weight","length","girth","pelage","days"))
cuweights=rbind(cuweights,cuweights.unmark,cuweights.dup)
#
# Define EN field
#
cuweights$EN=0
cuweights$EN[cuweights$cohort%in%ENYears]=1
if(any(is.na(cuweights$weight)))
{
	warning("\n Problem with data extraction. Following records are missing a weight. They have been removed.\n")
	for(w in which(is.na(cuweights$weight)))
		message("Animal ID = ",cuweights$AnimalID[w],"\n")
	cuweights=cuweights[!is.na(cuweights$weight),]
}

return(cuweights)
}

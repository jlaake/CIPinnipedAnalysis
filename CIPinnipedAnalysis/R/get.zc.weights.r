#' Extracts weights for Zc pups
#' :from various tables in the BrandMaster.mdb ACCESS database and returns a dataframe with ancillary fields added.
#' 
#' @param fdir directory for data files; if NULL uses location specified in databases.txt of CalcurData package; if "" uses databases in CalcurData pacakge; otherwise uses specified directory location
#' @param ENYears values of years that will be flagged as ENSO years
#' @param exclude if TRUE, it excludes any records with missing weight value; setting exclude=FALSE extracts all records which is useful for creating age-length key
#' @return dataframe of weights for Zc
#' @export
#' @author Jeff Laake
#' @examples 
#' # following runs adjustment model and all subsequent growth models
#' # and then stores table in CIPinnipedCensusQuery database. This will take
#' # a substantial amount of time
#' \donttest{
#' source(file.path(system.file(package="CIPinnipedAnalysis"),"ZcWeightAnalysis.r"))
#' }
get.zc.weights=function(fdir=NULL,ENYears=c(1976,1983,1984,1986,1987,1992,1997,1998,2002,2009),exclude=TRUE)
{
# read in area codes from Zc database
areas=getCalcurData("Zc","areacodes",dir=fdir)
#
# read in weights from ZcBrand table of the Access database
#
zcweights=getCalcurData("Zc","ZCBRAND",dir=fdir)
maxyear=max(zcweights$cohort)
zcweights$AnimalID=as.character(zcweights$AnimalID)
#
#  Exclude those with missing weight and unknown sex
#
if(exclude)
{
	zcweights=zcweights[!zcweights$sex=="U" & !is.na(zcweights$weight)& !is.na(zcweights$sex),]
} else
{
	zcweights=zcweights[!zcweights$sex=="U" & !is.na(zcweights$sex),]
}
zcweights$sex=factor(zcweights$sex)
#
#  read in the weights from the unmarkedpupweights table.  Exclude those with unknown sex or weight and from SCI and SNI.
#
zcweights.unmarked=getCalcurData("Zc","UnmarkedPupWeights",dir=fdir)
zcweights.unmarked=zcweights.unmarked[!zcweights.unmarked$sex=="U"& !is.na(zcweights.unmarked$sex),]
zcweights.unmarked$sex=factor(zcweights.unmarked$sex)
#
#  read in the tag only pup weights but exclude those with missing sex or weight
#
zcTagOnly=getCalcurData("Zc","TagInitial",dir=fdir)
if(exclude)
{
	zcTagOnly=zcTagOnly[zcTagOnly$age=="P"&!is.na(zcTagOnly$age)&!zcTagOnly$sex=="U" & !is.na(zcTagOnly$weight)& !is.na(zcTagOnly$sex),]
}else
{
	zcTagOnly=zcTagOnly[zcTagOnly$age=="P"&!is.na(zcTagOnly$age)&!zcTagOnly$sex=="U" & !is.na(zcTagOnly$sex),]
}
zcTagOnly$sex=factor(zcTagOnly$sex)
zcTagOnly=merge(zcTagOnly,subset(areas,select=c("region","sitecode")),by="region",all.x=TRUE)
#
#  read in the brand recaptures
#
zcRecap=getCalcurData("Zc","Recaptures",dir=fdir)
zcRecap$AnimalID=as.character(zcRecap$AnimalID)
if(exclude)
{
	zcRecap=zcRecap[!zcRecap$sex=="U" & !zcRecap$sex=="" &!is.na(zcRecap$weight)& !is.na(zcRecap$sex),]
}else
{
	zcRecap=zcRecap[!zcRecap$sex=="U" & !zcRecap$sex=="" & !is.na(zcRecap$sex),]
}
zcRecap$sex=factor(zcRecap$sex)
zcRecap=merge(zcRecap,subset(zcweights,select=c("cohort","AnimalID")),by="AnimalID")
#
#  read in the tag recaptures
#
zcTagRecap=getCalcurData("Zc","TagRecapture",dir=fdir)
if(exclude)
{
	zcTagRecap=zcTagRecap[!zcTagRecap$Sex=="U" & !is.na(zcTagRecap$Weight)& !is.na(zcTagRecap$Sex),]
}else
{
	zcTagRecap=zcTagRecap[!zcTagRecap$Sex=="U" & !is.na(zcTagRecap$Sex),]
}
zcTagRecap$sex=factor(zcTagRecap$Sex)
zcTagRecap=merge(zcTagRecap,subset(zcTagOnly,select=c("cohort","AnimalID")),by="AnimalID")
#
#   Add a days field which is the number of days from 1 Oct of the year
#
zcweights$days=floor(as.numeric((zcweights$branddate-as.POSIXct(paste(as.character(zcweights$cohort),"-10-01",sep=""),format="%Y-%m-%d"))/(24*3600)))
zcTagOnly$days=floor(as.numeric((zcTagOnly$capturedate-as.POSIXct(paste(as.character(zcTagOnly$cohort),"-10-01",sep=""),format="%Y-%m-%d"))))
zcweights.unmarked$days=floor(as.numeric((zcweights.unmarked$capturedate-as.POSIXct(paste(as.character(zcweights.unmarked$cohort),"-10-01",sep=""),format="%Y-%m-%d"))/(24*3600)))
zcRecap$days=floor(as.numeric((zcRecap$capturedate-as.POSIXct(paste(as.character(zcRecap$cohort),"-10-01",sep=""),format="%Y-%m-%d"))))
zcTagRecap$days=floor(as.numeric((zcTagRecap$capturedate-as.POSIXct(paste(as.character(zcTagRecap$cohort),"-10-01",sep=""),format="%Y-%m-%d"))))
zcRecap=zcRecap[zcRecap$days>0,]
#
#  Combine data tables
#
zcweights.unmarked=subset(zcweights.unmarked,select=c("AnimalID","sitecode","cohort","sex","weight","length","girth","days"))
if(exclude) zcweights.unmarked=zcweights.unmarked[!is.na(zcweights.unmarked$weight),]
zcweights$sitecode="SMI"
zcweights=subset(zcweights,select=c("AnimalID","sitecode","cohort","sex","weight","length","girth","days"))
#zcTagOnly$sitecode="SMI"
zcTagOnly=subset(zcTagOnly,select=c("AnimalID","sitecode","cohort","sex","weight","length","girth","days"))
zcRecap=subset(zcRecap,select=c("AnimalID","sitecode","cohort","sex","weight","length","girth","days"))
#zcTagRecap$sitecode="SMI"
zcTagRecap=subset(zcTagRecap,select=c("AnimalID","sitecode","cohort","sex","Weight","Length","Girth","days"))
names(zcTagRecap)=c("AnimalID","sitecode","cohort","sex","weight","length","girth","days")
zcweights$AnimalID=as.character(zcweights$AnimalID)
zcTagOnly$AnimalID=as.character(zcTagOnly$AnimalID)
zcTagRecap$AnimalID=as.character(zcTagRecap$AnimalID)
zcRecap$AnimalID=as.character(zcRecap$AnimalID)
zcweights.unmarked$AnimalID=as.character(zcweights.unmarked$AnimalID)
zcweights$tag.type="Brand"
zcweights.unmarked$tag.type="Unmarked"
zcRecap$tag.type="Brand"
zcTagRecap$tag.type="Tag"
zcTagOnly$tag.type="Tag"
zcweights$type="Initial"
zcweights.unmarked$type="Initial"
zcTagOnly$type="Initial"
zcTagRecap$type="Recap"
zcRecap$type="Recap"
zcweights=rbind(zcTagOnly,zcweights,zcweights.unmarked,zcRecap,zcTagRecap)
zcweights$type=factor(zcweights$type)
zcweights$sitecode=factor(zcweights$sitecode)
zcweights$tag.type=factor(zcweights$tag.type)
#
# Define EN field
#
zcweights$EN=0
zcweights$EN[zcweights$cohort%in%ENYears]=1
#
# check for any NA weights
#
if(exclude & any(is.na(zcweights$weight)))
{
	warning("\n Problem with data extraction. Following records are missing a weight. They have been removed.\n")
	for(w in which(is.na(zcweights$weight)))
		message("Animal ID = ",zcweights$AnimalID[w],"\n")
	zcweights=zcweights[!is.na(zcweights$weight),]
}
return(zcweights)
}





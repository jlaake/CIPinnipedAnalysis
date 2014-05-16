#' Checks for errors in the assigned dead pup sample areas and CU Survey dates
#' 
#' Checks to make sure that each table in CIPinnipedCensusMaster.mdb has a
#' matching record in DeadPupSampleAreas table for SMI for species=Zc and Cu.
#' 
#' Also reports any CU survey dates in dead file that are more than 5 days from
#' assigned CU Survey Date for that survey number in that year.  Typically the error is 
#' 
#' \preformatted{ If there are no missing records the output should look as
#' follows:
#' ***Checking Zc CU dead pup census for Zc at SMI ***Checking Zc CU dead pup
#' census for Cu at SMI ***Checking Zc CU live pup census for Zc at SMI
#' ***Checking Zc CU live pup census for Cu at SMI ***Checking Zc dead tag
#' initial for Zc at SMI ***Checking Zc dead tag resight for Zc at SMI } If
#' there are records that don't have a matching record in DeadPupSampleAreas
#' they are printed out. DeadPupSampleAreas must contain a record for each Area
#' Code for each year that the Area Code was used in either the tables for dead
#' pups or live pups because that tells the code how to sum the live and dead
#' counts for computation of the mortality rate.
#' 
#' @export
#' @import RODBC
#' @param fdir directory for CIPinnipedCensusMaster.mdb file
#' @return None
#' @note Creates log file MissingDeadPupArea.txt with any errors
#' @author Jeff Laake 
check.dead.pup.areas <-
function(fdir=NULL)
{
#
#  Make connection to CIPinnipedCensusMaster.mdb
if(!is.null(fdir) && fdir=="")
{
	fdir=system.file(package="CalcurData")
} 
#
sink("MissingDeadPupArea.txt")
areas=getCalcurData("CIPCensus","DeadPupSampleAreas",dir=fdir)
areas$YearArea=as.character(areas$YearArea)
dead=getCalcurData("CIPCensus","Zc Cu dead pup census",dir=fdir)
dead$key=paste(dead$Year,dead[,"Survey number"])
cusurvey=getCalcurData("CIPCensus","CU Survey Dates",dir=fdir)
cusurvey$key=paste(cusurvey$Year,cusurvey[,"Survey number"])
live=getCalcurData("CIPCensus","Zc Cu live pup census",dir=fdir)
taginitial=getCalcurData("CIPCensus","Zc dead tag initial",dir=fdir)
tagresight=getCalcurData("CIPCensus","Zc dead tag resight",dir=fdir)
dead_withdates=merge(dead[dead$Species=="Cu"&!dead[,"Area code"]%in%c("WCR","ECR"),],cusurvey,by="key",all.x=TRUE)
cat("\n\n ***Checking CU Survey Dates against CU survey dates for SMI mainland\n")
if(any(is.na(dead_withdates$Date)))
{
  cat("\nMissing survey numbers in CU Survey Dates for the following: \n")
  print(dead_withdates[is.na(dead_withdates$Date),])  
}
dead_withdates=dead_withdates[!is.na(dead_withdates[,"Survey date"]),]
date_gaps=(dead_withdates[,"Survey date"]-dead_withdates$Date)/(24*3600)
if(any(abs(date_gaps)>5))
{
	cat("\nDate gap between survey date in dead file and date in CU Survey > 5 days \n")
	dead_withdates=dead_withdates[abs(date_gaps)>5,c(4,6,21)]
	names(dead_withdates)[1]=c("Survey date in dead file")
	names(dead_withdates)[3]=c("Survey date in CU Survey Date file")
	print(dead_withdates) 
}
cat("\n\n ***Checking Zc CU dead pup census for Zc at SMI\n")
x=dead[dead$Species=="Zc"&dead$Island=="SMI",]
x$YearArea=paste(x$Year,x[,"Area code"],sep="")
z=merge(x,areas[areas$YearAreaSpecies=="Zc",],all.x=TRUE,by="YearArea")
if(any(is.na(z[,"Dead pup sample area"])))
{
   cat("\nMissing DeadPupSampleAreas for Zc at SMI\n")
   print(z[is.na(z[,"Dead pup sample area"]),])
}
cat("\n\n ***Checking Zc CU dead pup census for Cu at SMI\n")
x=dead[dead$Species=="Cu"&dead$Island=="SMI",]
x$YearArea=paste(x$Year,x[,"Area code"],sep="")
z=merge(x,areas[areas$YearAreaSpecies=="Cu",],all.x=TRUE,by="YearArea")
if(any(is.na(z[,"Dead pup sample area"])))
{
   cat("\nMissing DeadPupSampleAreas for Cu at SMI\n")
   print(z[is.na(z[,"Dead pup sample area"]),])
}
cat("\n\n ***Checking Zc CU live pup census for Zc at SMI\n")
x=live[live$Species=="Zc",]
x$YearArea=paste(x$Year,x[,"Area code"],sep="")
z=merge(x,areas[areas$YearAreaSpecies=="Zc",],all.x=TRUE,by="YearArea")
if(any(is.na(z[,"Dead pup sample area"])))
{
   cat("\nMissing DeadPupSampleAreas for Zc at SMI\n")
   print(z[is.na(z[,"Dead pup sample area"]),])
}
cat("\n\n ***Checking Zc CU live pup census for Cu at SMI\n")
x=live[live$Species=="Cu",]
x$YearArea=paste(x$Year,x[,"Area code"],sep="")
z=merge(x,areas[areas$YearAreaSpecies=="Cu",],all.x=TRUE,by="YearArea")
if(any(is.na(z[,"Dead pup sample area"])))
{
   cat("\nMissing DeadPupSampleAreas for Cu at SMI\n")
   print(z[is.na(z[,"Dead pup sample area"]),])
}
cat("\n\n ***Checking Zc dead tag initial for Zc at SMI\n")
x=taginitial[taginitial$Island=="SMI",]
x$YearArea=paste(x$Year,x[,"Area code"],sep="")
z=merge(x,areas[areas$YearAreaSpecies=="Zc",],all.x=TRUE,by="YearArea")
if(any(is.na(z[,"Dead pup sample area"])))
{
   cat("\nMissing DeadPupSampleAreas for Zc at SMI\n")
   print(z[is.na(z[,"Dead pup sample area"]),])
}
cat("\n\n ***Checking Zc dead tag resight for Zc at SMI\n")
x=tagresight[tagresight$Island=="SMI",]
x$YearArea=paste(x$Year,x[,"Area code"],sep="")
z=merge(x,areas[areas$YearAreaSpecies=="Zc",],all.x=TRUE,by="YearArea")
if(any(is.na(z[,"Dead pup sample area"])))
{
   cat("\nMissing DeadPupSampleAreas for Zc at SMI\n")
   print(z[is.na(z[,"Dead pup sample area"]),])
}
sink()
invisible()
}


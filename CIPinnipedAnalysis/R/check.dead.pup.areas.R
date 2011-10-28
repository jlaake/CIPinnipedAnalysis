

#' Checks for errors in the assigned dead pup sample areas
#' Checks to make sure that each table in CIPinnipedCensusMaster.mdb has a
#' matching record in DeadPupSampleAreas table for SMI for species=Zc and Cu.
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
function(fdir="")
{
#
#  Make connection to CIPinnipedCensusMaster.mdb
#
if(fdir=="")fdir=system.file(package="CIPinnipedAnalysis")
fdir=file.path(fdir,"Master/CIPinnipedCensusMaster.mdb")
connection=odbcConnectAccess2007(fdir)
sink("MissingDeadPupArea.txt")
areas=sqlFetch(connection,"DeadPupSampleAreas")
areas$YearArea=as.character(areas$YearArea)
dead=sqlFetch(connection,"Zc Cu dead pup census")
live=sqlFetch(connection,"Zc Cu live pup census")
taginitial=sqlFetch(connection,"Zc dead tag initial")
tagresight=sqlFetch(connection,"Zc dead tag resight")
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
odbcCloseAll()
sink()
invisible()
}


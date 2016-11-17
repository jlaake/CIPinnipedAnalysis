#
#  Constructs Pup Production table for  Zc at SMI and SNI in CIPinnipedCensusQuery
#  using live and dead pup counts from CIPinnipedCensusMaster via links and queries
#  in CIPinnipedCensusQuery.  
# 
# Uses the production.stats to create ZcProduction table for SMI and SNI in CIPinnipedCensusQuery.mdb.
# Any expected warnings or error messages are printed in PupProduction.log and unexpected
# errors/coding problems would be found in PupProduction.out.  If there are
# any errors then the production tables will not be created. 
# 
# fdir directory for data files; if NULL uses location specified in databases.txt of CalcurData package; if "" uses databases in CalcurData pacakge; otherwise uses specified directory location
if(!exists("fdir"))fdir=NULL
sink("PupProduction.log")
if(!is.null(fdir) && fdir=="")
{
	fdir1=system.file(package="CalcurData")
	fdir2=fdir1
} else
{
	fdir1=fdir
	fdir2=fdir
}
# Compute production stats for smi and castle rock for Zc
smidat=production.stats(island="SMI",mainland=TRUE,species="Zc",dir=fdir2)
crdat=production.stats(island="SMI",mainland=FALSE,species="Zc",dir=fdir2)
# Compute production stats for San Nicolas
livepups=getCalcurData("CIPCensus","Zc Cu live pup census",dir=fdir2)
snidat=production.stats(island="SNI",mainland=FALSE,species="Zc",dir=fdir2,years=sort(unique(livepups$Year[livepups$Island=="SNI"])))
rm(livepups)
xx=saveCalcurData(rbind(smidat,crdat,snidat),db="CIPquery",tbl="ZcProduction",dir=fdir1)
sink()



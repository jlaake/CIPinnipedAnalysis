#
#  Constructs Pup Production tables for Cu and Zc at SMI in CIPinnipedCensusQuery
#  using live and dead pup counts from CIPinnipedCensusMaster via links and queries
#  in CIPinnipedCensusQuery.  It also saves plots of production in pdf files.
#
# 
# Uses the production.stats to create ZcProduction and CuProduction tables for SMI in CIPinnipedCensusQuery.mdb.
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
	sdir=system.file(package="CIPinnipedAnalysis")
	source(file.path(sdir,"ZcPupProduction.r"))
	source(file.path(sdir,"CuPupProduction.r"))	
	sink()

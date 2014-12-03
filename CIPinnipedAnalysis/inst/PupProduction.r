#
#  Constructs Pup Production tables for Cu and Zc at SMI in CIPinnipedCensusQuery
#  using live and dead pup counts from CIPinnipedCensusMaster via links and queries
#  in CIPinnipedCensusQuery.  It also saves plots of production in pdf files.
#
# 
# For Zc and Cu, creates data tables in the ACCESS CIPinnipedCensusQuery database with the early pup
# mortality estimates during the season in each year.  Also constructs pdf
# plots of those values.
# 
# Uses the production.stats to create ZcProduction and CuProduction tables for SMI in CIPinnipedCensusQuery.mdb.
# Any expected warnings or error messages are printed in PupProduction.log and unexpected
# errors/coding problems would be found in PupProduction.out.  If there are
# any errors then the production tables will not be created. In addition to
# creating the tables, it also produces a set of plots of the results in pdfs.
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
# Compute production stats for CU on SMI and Castle Rock
	smidat=production.stats(island="SMI",mainland=TRUE,species="Cu",dir=fdir2)
	crdat=production.stats(island="SMI",mainland=FALSE,species="Cu",dir=fdir2)
	xx=saveCalcurData(rbind(smidat,crdat),db="CIPquery",tbl="CuProduction",dir=fdir1)
	sink()
# Get production tables and create pdfs
	ZcProduction=getCalcurData("CIPquery","ZcProduction",dir=fdir1)
	CuProduction=getCalcurData("CIPquery","CuProduction",dir=fdir1)


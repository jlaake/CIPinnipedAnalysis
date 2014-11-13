#' Store weights in CIPinnipedCensusQuery
#' 
#' Stores standardized weights into table in ACCESS database for Zc and Cu
#' depending on the value of species.
#' 
#' @param x dataframe containing weights to be stored in table
#' @param species species code; either Zc or Cu
#' @param fdir directory for databases; if NULL uses location from databases.txt in CalcurData package; if "" uses databases stored in package directory; if neither you can specify the location
#' @return None
#' @export
#' @author Jeff Laake
#' 
store_weights=function(x,species="Zc",fdir=NULL)
{
    # Export to CIPinnipedCensusQuery
	years=as.numeric(rownames(x))
	# this code takes table of weights which are recorded with lots of digits and 
	# converts to 3 digits of precision; the warnings are suppressed because some values
	# are missing and NA will generate a warning
	suppressWarnings(
			xx<-apply(x[,-1],2,function(x) as.numeric(sprintf("%.3f", x))))
	xf=data.frame(Year=years,cbind(as.data.frame(xx)))
	# if Zc save in ZcWeights otherwise in CuWeights
	if(species=="Zc")
	   xx=saveCalcurData(x,db="CIPquery",tbl="ZcWeights",dir=fdir)
	else
	xx=saveCalcurData(CUWeight.df,db="CIPquery",tbl="CuWeights",dir=fdir)
    return(NULL)	
}

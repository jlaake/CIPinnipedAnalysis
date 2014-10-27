#' Store weights in CIPinnipedCensusQuery
#' 
#' Stores standardized weights into table in ACCESS database for Zc and Cu
#' depending on the value of species.
#' 
#' @param x dataframe containing weights to be stored in table
#' @param species species code; either Zc or Cu
#' @return None
#' @export
#' @author Jeff Laake
#' 
store_weights=function(x,species="Zc")
{
    # Export to CIPinnipedCensusQuery
	years=as.numeric(rownames(x))
	suppressWarnings(
			xx<-apply(x[,-1],2,function(x) as.numeric(sprintf("%.3f", x))))
	xf=data.frame(Year=years,cbind(as.data.frame(xx)))
	if(species=="Zc")
	   xx=saveCalcurData(x,db="CIPquery",tbl="ZcWeights",dir=fdir)
	else
	xx=saveCalcurData(CUWeight.df,db="CIPquery",tbl="CuWeights",dir=fdir)
    return(NULL)	
}

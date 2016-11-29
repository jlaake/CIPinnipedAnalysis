#' Create average anomalies for a range of months 
#' 
#' @param v matrix of environmental values with columns for months and rows for years
#' @param start first location to use when v is treated as a vector of the matrix transposed
#' @param n number of months to include in each average
#' @return vector of averages with year names
#' @export
#' @author Jeff Laake
average_anomalies=function(v,start,n)
{
	ir=ceiling(start/12)
	years=rownames(v)[ir:nrow(v)]
	v=as.vector(t(v))
	if(start>1) v=v[-(1:(start-1))]
	avg=apply(cbind(seq(1,length(v),12),(seq(1,length(v),12)+n-1)),1,function(x)mean(v[x[1]:x[2]],na.rm=TRUE))
	names(avg)=years
	return(avg)
}

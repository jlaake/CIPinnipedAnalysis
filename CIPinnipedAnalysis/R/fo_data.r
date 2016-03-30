#' Creates Frequency of Occurrence (FO) DataFrame for Primary Prey
#' 
#' @return dataframe of FO for key prey or prey groupings and derived values
#' including first and second principal components of 5 primary prey
#' @export
#' @author Jeff Laake
#'
create_fo=function()
{
	data(fo)
	sardine=get_fo(fo,c("SARSAG"))$fo
	anchovy=get_fo(fo,c("ENGMOR"))$fo
	sa=get_fo(fo,c("SARSAG","ENGMOR"))$fo
	rockfish=get_fo(fo,c("SEBSPP","SEBSP1","SEBSP2"))$fo
	hake=get_fo(fo,c("MERPRO"))$fo
	squid=get_fo(fo,c("LOLOPA"))$fo
	mackerel=get_fo(fo,c("SCOJAP","SCOSPP","SCOMAR","TRASYM"))$fo
	lo.cal=get_fo(fo,c("SEBSPP","SEBSP1","SEBSP2","MERPRO","LOLOPA"))$fo
	hi.cal=get_fo(fo,c("SARSAG","ENGMOR","SCOJAP","SliCOSPP","SCOMAR","TRASYM"))	
	fo=cbind(Year=as.numeric(names(sardine)),sardine=sardine,anchovy=anchovy,
			rockfish=rockfish,hake=hake,squid=squid,mackerel=mackerel)
# first and second principal components of 5 primary prey species
	fo=cbind(fo,predict(prcomp(fo[,-1]))[,1:2])
	fo=as.data.frame(fo)
	fo=cbind(fo,sa=as.vector(sa),lo.cal=as.vector(lo.cal),hi.cal=as.vector(hi.cal$fo),n.scat=as.vector(hi.cal$n.scat))
	fo$ratio=fo$hi.cal/fo$lo.cal
	fo$diet=NA
	fo$diet[fo$Year%in%1981:1986]="Diet1"
	fo$diet[fo$Year%in%c(1980,2009,2012,2013,1992,2000,1995,2001,2011,1991,2010)]="Diet2"
	fo$diet[fo$Year%in%c(1993,1996,2002,2003,2004,2005,2006,2007)]="Diet3"
	fo$diet=factor(fo$diet)
	return(fo)
}
#' Creates Frequency of Occurrence (FO) DataFrame for Primary Prey of Callorhinus
#' 
#' @return dataframe of FO for key prey or prey groupings 
#' @export
#' @author Jeff Laake
#'
create_cufo=function(fo)
{
	scats=NULL
	sardine=get_fo(fo,c("SARSAG"))$fo
	anchovy=get_fo(fo,c("ENGMOR","ENGSPP"))$fo
	rockfish=get_fo(fo,c("SEBSPP","SEBSP1","SEBSP2"))$fo
	hake=get_fo(fo,c("MERPRO"))$fo
	squid=get_fo(fo,c("LOLOPA"))$fo
	mackerel=get_fo(fo,c("TRASYM"))$fo
	octupus=get_fo(fo,c("OCTBIM","OCTRUB"))$fo
	leusti=get_fo(fo,c("LEUSTI"))$fo
	cephalopod=get_fo(fo,c("ABRFEL","ONYBOR","GONSPP","GONSPB","GONONY","GONBOR","GONPYR","GONSPX"))$fo
	myctophid=get_fo(fo,c("SYMCAL"))$fo
	fo=cbind(Year=as.numeric(names(sardine)),sardine=sardine,anchovy=anchovy,
			rockfish=rockfish,hake=hake,squid=squid,mackerel=mackerel,octupus=octupus,cephalopod=cephalopod,
			leusti=leusti,myctophid=myctophid,n.scat=sapply(split(fo$SCATNUM,fo$Year),function(x) length(unique(x))))
	return(fo)
}
#' Computes Frequency of Occurrence (FO) for Diet Data
#' 
#' For a given set of species codes it computes FO by extracting the observed scats and
#' combines into year groupings. For years before 1992 when only otoliths were
#' used, it "corrects" based on data collected after 1991. 
#' 
#' Dataframe scats no longer used because some had entries that were not identified so fo was too large.  Now the counts
#' of scats are calculated from unique number of scats in prey (fo). This means for the 1980s any scat that could have 
#' had bone only is no longer in the count, so the bone_cf is not used or fo values greater than one could be obtained.
#' Thus up through 1991, fo is based on otloliths and beaks only.
#' 
#' @param fo dataframe of prey species identified in each scat and whether it was identified by otolith/beak or other 
#' @param species a vector of prey species codes for which FO is computed. For more than one species code
#' it is the FO of any of the species.
#' @param exclude if TRUE, then uses any species other than those listed in species
#' @param years a vector of cut points for year groupings
#' @param labels to be used for year groupings
#' @return A list with vector of FO values for each year in the data and number of scats for each year
#' @export
#' @author Jeff Laake
#'
get_fo=function(fo,species,exclude=FALSE,years=NULL,labels=NULL)
{
	if(exclude)
		fosubset=fo[!fo$PREYSP%in%species,]
	else
		fosubset=fo[fo$PREYSP%in%species,]
	cf=bone_cf(fosubset)
	
	if(!is.null(years))
		fo$Year=cut(fo$Year,years,labels=labels)
	else
		fo$Year=factor(fo$Year)
	size=sapply(split(fo$SCATNUM,fo$Year),function(x) length(unique(x)))
	if(is.null(years))
		fosubset$Yr=factor(fosubset$Year,levels=names(size))	
	else
		fosubset$Yr=cut(fosubset$Year,years,labels=labels)
	#fosubset$fo=ifelse(fosubset$Year<=1991,cf$cf,1)
	fosubset$fo=1
	fotable=tapply(fosubset$fo,list(fosubset$Yr,fosubset$SCATNUM),mean)
	fotable[is.na(fotable)]=0
	return(list(fo=rowSums(fotable)/size,n.scat=size))
}
#' Computes correction factor for pre-1992 data
#' 
#' For a given set of species codes it computes the correction factor for
#' oto/bk only data collected before 1992. It computes FO with all structures divided by
#' FO from oto/bk only and that ratio is then used to correct years prior to 1992 in
#' the fo function. 
#' 
#' No longer used because scats are counted as the number of scats with any identified prey from the 
#' fo dataframe.
#' 
#' @param fo dataframe of prey species identified in each scat and whether it was identified by otolith/beak or other 
#' @return list with correction factor cf and its std error se.cf
#' @export
#' @author Jeff Laake
#'
bone_cf=function(fo)
{
	fo$Year=factor(fo$Year)
	X=apply(with(fo[fo$FO_OT.BK==1,],table(Year,SCATNUM)),1,function(x) sum(x>0))
	Y=apply(with(fo,table(Year,SCATNUM)),1,function(x) sum(x>0))
	Y=Y[as.numeric(names(Y))>=1992]
	X=X[as.numeric(names(X))>=1992]
	meanX=mean(X)
	count=length(X)
	cf=sum(Y)/sum(X)
	vary=var(Y)
	varx=var(X)
	covxy=cov(X,Y)
	se.cf=sqrt((vary+cf^2*varx-2*cf*covxy)/(count*meanX^2))
	return(list(cf=cf,se.cf=se.cf))
}

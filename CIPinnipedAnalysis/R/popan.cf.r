#' Construct correction factor data from dead pup tagging studies
#' 
#' For a given island (SMI or SNI) and year, the dead pup tagging data is analyzed with a POPAN 
#' model to create correction factors for the observed dead pup counts to estimate the total number that died.
#' 
#' The function getdead_ch is used to extract the data from the CIPinnipedCensusMaster database and 
#' create capture histories for tagged and untagged pups. This is then analyzed with the POPAN model in MARK using
#' the RMark package.  A set of models is examined that depends on the year because the number of sampling occasions
#' and the data collected has varied.  Prior to 1998, the position on the beach and substrate were not recorded for
#' dead pups when they were stacked.  However, the area was recorded and the substrate is largely similar in each area. Thuse
#' for years prior to 1998 the area was used as a grouping factor and in 1998 and beyond, substrate and position were used as
#' grouping factors. In POPAN models, p on the first and last occasion are not separately estimable in a time model
#' so these have been restricted such that p for occasion 1 is the same as occasion 2 and p for the last occasion is the same as p for the 
#' next to last occasion (see definition of time2). 
#' 
#' Each set of models is run and then they are re-run using the initial values from the best model.  This helps with
#' model convergence which can fail for POPAN models.  To see which models are failing set silent=FALSE. Once the models have been
#' run, the function popan.derived is called with the list of model results which computes the abundances and related statistics (immigration BiGross) by model averaging
#' over the set of models.  These are then stored in the list cfdata which is returned with the correction factors by occasion (cfbyocc) and the marklist of final model results.
#' @export  
#' @param island ("SMI" or "SNI")
#' @param year four digit numeric year
#' @param silent if TRUE shows each model as it is run and any problems that occur; if FALSE this is hidden each survey date (occasion); if non-null provides correction at that date
#' @return a list containing a marlist of the final model results, a list with the correction factor data (cfdata) and correction factors by occasion (cfbyocc). 
#' @author Jeff Laake 
popan.cf=function(island,year,silent=TRUE)
{
	
	x.popan=getdead_ch(island,year)
	
	chlist=CIPinnipedAnalysis:::process_ch(x.popan$df$ch,x.popan$df$freq)
	
	times=x.popan$days
	time.intervals=as.numeric(diff(times))
	if(year>1997)
		x.proc=process.data(x.popan$df,model="POPAN",time.intervals=time.intervals,groups=c("Position","Substrate"),begin.time=0)
	else
		x.proc=process.data(x.popan$df,model="POPAN",time.intervals=time.intervals,groups=c("Area code"),begin.time=0)
	
	numdead= tapply(abs(x.proc$data$freq),list(factor(chlist$first,1:nrow(times)),x.proc$data$group),sum)
	numdead[is.na(numdead)]=0
	
	x.ddl=make.design.data(x.proc)
	x.ddl$p$time2=x.ddl$p$time
	x.ddl$p$time2[x.ddl$p$time==times[2]]=times[1]
	if(year%in%c(1998,2002))
	{
		x.ddl$p$time2[x.ddl$p$time==times[3]]=times[1]
		x.ddl$p$time2[x.ddl$p$time==times[4]]=times[1]
	}
	x.ddl$p$time2[x.ddl$p$time==times[length(times)-1]]=times[length(times)]
	x.ddl$p$time2=factor(as.character(x.ddl$p$time2))
	
	x.ddl$Phi$time2=x.ddl$Phi$time
	if(year%in%c(1998,2002))
	{
		x.ddl$Phi$time2[x.ddl$Phi$time==times[2]]=times[1]
		x.ddl$Phi$time2[x.ddl$Phi$time==times[3]]=times[1]
		x.ddl$Phi$time2=factor(as.character(x.ddl$Phi$time2))
	}
	
	
	if(year>1997)
	{
		do.cf=function(initial=NULL)
		{
			Phi.1=list(formula=~Substrate+Position)
			Phi.2=list(formula=~Substrate)
			Phi.3=list(formula=~Position)
			Phi.4=list(formula=~1)
			Phi.5=list(formula=~Substrate+Position+time2)
			Phi.6=list(formula=~Substrate+time2)
			Phi.7=list(formula=~Position+time2)
			Phi.8=list(formula=~time2)
			Phi.9=list(formula=~Substrate+Position+Time)
			Phi.10=list(formula=~Substrate+Time)
			Phi.11=list(formula=~Position+Time)
			Phi.12=list(formula=~Time)
			pent.1=list(formula=~time)
			pent.2=list(formula=~1)
			pent.3=list(formula=~Time)
			p.1=list(formula=~time2+Substrate)
			p.2=list(formula=~Substrate)
			p.3=list(formula=~time2)
			p.4=list(formula=~1)
			N.1=list(formula=~group)
			cml=create.model.list("POPAN")
			results=mark.wrapper(cml,data=x.proc,ddl=x.ddl,initial=initial,output=FALSE,silent=silent)
			return(results)
		}
	}else
	if(year<1997)
	{
		do.cf=function(initial=NULL)
		{
			Phi.1=list(formula=~group)
			Phi.2=list(formula=~1)
			Phi.3=list(formula=~group+time2)
			Phi.4=list(formula=~time2)
			Phi.5=list(formula=~group+Time)
			Phi.6=list(formula=~Time)
			pent.1=list(formula=~time)
			pent.2=list(formula=~1)
			pent.3=list(formula=~Time)
			p.1=list(formula=~time2+group)
			p.2=list(formula=~group)
			p.3=list(formula=~time2)
			p.4=list(formula=~1)
			N.1=list(formula=~group)
			cml=create.model.list("POPAN")
			results=mark.wrapper(cml,data=x.proc,ddl=x.ddl,initial=initial,output=FALSE,silent=silent)
			return(results)
		}
	}
	cat("\nrunning popan models\n")
	popan.results=do.cf()
	cat("\nre-running models with initial values\n")
	popan.results=do.cf(popan.results[[as.numeric(row.names(popan.results$model.table[1,]))]])
		
	cfdata=suppressMessages(popan.derived(x.proc,popan.results))
		
	cfdata$BiGross$daysfrom1July=rep(x.popan$daysfrom1July,nrow(x.proc$group.covariates))
	cfdata$BiGross$cumdead=as.vector(apply(t(numdead), 1, cumsum))
	cfdata$BiGross$cf=max(1,cfdata$BiGross$estimate/cfdata$BiGross$cumdead)
	cfdata$BiGross$cf.se=cfdata$BiGross$se/cfdata$BiGross$cumdead
	cfbyocc=with(cfdata$BiGross,tapply(estimate,occasion,sum)/tapply(cumdead,occasion,sum))
	cfbyocc=data.frame(occasion=names(cfbyocc),daysfrom1July=x.popan$daysfrom1July,cfbyocc)
		
	k=length(times)
	vcmat=cbind(1:nrow(x.proc$group.covariates),rep(1:k,each=nrow(cfdata$BiGross.vcv)),as.vector(cfdata$BiGross.vcv))
	se=rep(0,k)
	for(i in 1:k)
	{
		se[i]=sqrt(sum(vcmat[vcmat[,1]==i | vcmat[,2]==i,]))
	} 
	cfbyocc$se=se/tapply(cfdata$BiGross$cumdead,cfdata$BiGross$occasion,sum)
	cfbyocc$observed.dead=tapply(cfdata$BiGross$cumdead,cfdata$BiGross$occasion,sum)
	cfbyocc$estimated.dead=with(cfdata$BiGross,tapply(estimate,occasion,sum))
		
	rownames(cfdata$NGross)= apply(x.proc$group.covariates, 1, paste, collapse = "")
		
	return(list(popan.results=popan.results,cfdata=cfdata,cfbyocc=cfbyocc))
	}
	
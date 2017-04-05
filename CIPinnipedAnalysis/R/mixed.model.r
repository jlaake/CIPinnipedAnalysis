#' Fit a sequence of mixed effect models with normal error distribution
#' 
#' Either fits a sequence of random effect models with a single fixed effect formula with REML or
#' fits a sequence of fixed effect models with a single random effect using MLE.
#' 
#' @param fixed.f sequence of formulas for fixed effects 
#' @param random.f sequence of formulas for random effects 
#' @param data dataframe for model fitting
#' @param control argument for nlme
#' @param method REML or ML; chosen as needed but can be specified
#' @param save.model if TRUE saves each complete model rather than removing data,fitted,residual components
#' @param ... other arguments for plot
#' @return list with following elements: results - list with each fitted model object; model.table- model selection table;
#' best - model number that was best fit in results list; best.f - formula for fixed effect from best model; best.r - formula for random effects for best model
#' @export
#' @author Jeff Laake
fitmixed=function(fixed.f,random.f,data,control=lmeControl(opt="optim"),method=NULL,save.model=FALSE,...)
{
	if(class(fixed.f)=="formula") fixed.f=list(fixed.f)
	if(length(fixed.f)>1)
	{
		if(length(random.f)>1) stop("cannot specify multiple sets of fixed and random formula")
		if(is.null(method))method="ML"
	} else
	{
		if(length(random.f)>1)
			if(is.null(method))method="REML"
	}
	# subset data to only include all variables used in the set of models and then exclude any rows with NA
	data=subset(data,select=unique(c(unlist(sapply(fixed.f,all.vars)),unlist(sapply(unlist(random.f),all.vars)))))
	data=na.exclude(data)
	results=vector("list",length(fixed.f)*length(random.f))
	model.table=data.frame(fixed=rep(NA,length(fixed.f)*length(random.f)),random=rep(NA,length(fixed.f)*length(random.f)))
	i=0
	for(f in fixed.f)
		for(r in random.f)
		{
			i=i+1
			results[[i]]=lme(f,random=r,data=data,control=control,method=method,...)
			if(!save.model)
			{
				results[[i]]$data=NULL
			    results[[i]]$fitted=NULL
			    results[[i]]$residuals=NULL
			}
			cf=as.character(f)
			cf[1]=cf[2]
			cf[2]="~"
			model.table[i,]=cbind(fixed=paste(cf,collapse=""),random=paste(r,collapse=", "))
		}
	model.table$AIC=sapply(results,AIC)
	results$model.table=model.table[order(model.table$AIC),]
	results$best=as.numeric(row.names(results$model.table)[1]) 
	if(length(fixed.f)==1)
		results$best.f=fixed.f[[1]]
	else
		results$best.f=fixed.f[[results$best]]
	if(length(random.f)==1)
		results$best.r=random.f[[1]]
	else
		results$best.r=random.f[[results$best]]
	results$data=data
	return(results)
}
#' Compute AICc (AIC with small sample size correction) for a series of models and model.table with
#' number of parameters and model weight.
#' 
#' @param results list of model results with element model.table
#' @param nobs number of observations  
#' @param nref number of random effect parameters
#' @param fixed.f vector of fixed effect formulas
#' @return list with following elements: results - list with each fitted model object; model.table- model selection table;
#' best - model number that was best fit in results list; best.f - formula for fixed effect from best model; best.r - formula for random effects for best model
#' @export
#' @author Jeff Laake

compute_AICc=function(results,nobs,nref,fixed.f)
{
	i=1
	results$model.table$AICc=NA
	results$model.table$fixed_par=NA
	results$model.table$random_par=NA
	for( mod in as.numeric(row.names(results$model.table)))
	{
		K=ncol(coef(results[[mod]]))+nref
		results$model.table$fixed_par[i]=ncol(coef(results[[mod]]))
		results$model.table$random_par[i]=nref
		results$model.table$AICc[i]=results$model.table$AIC[i] + 2*K*(K+1)/(nobs-K+1)
		i=i+1
	}
	results$model.table=results$model.table[order(results$model.table$AICc),]
	results$model.table$weight=exp(-.5*(results$model.table$AICc-min(results$model.table$AICc)))
	results$model.table$weight=results$model.table$weight/sum(results$model.table$weight)
	results$best=as.numeric(row.names(results$model.table))[1]
	results$best.f=fixed.f[[results$best]]
	return(results)
}

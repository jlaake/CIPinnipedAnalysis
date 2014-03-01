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
#' @param ... other arguments for plot
#' @return list with following elements: results - list with each fitted model object; model.table- model selection table;
#' best - model number that was best fit in results list; best.f - formula for fixed effect from best model; best.r - formula for random effects for best model
#' @export
#' @author Jeff Laake
fitmixed=function(fixed.f,random.f,data,control=lmeControl(opt="optim"),method=NULL,...)
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
	results=vector("list",length(fixed.f)*length(random.f))
	model.table=data.frame(fixed=rep(NA,length(fixed.f)*length(random.f)),random=rep(NA,length(fixed.f)*length(random.f)))
	i=0
	for(f in fixed.f)
		for(r in random.f)
		{
			i=i+1
			results[[i]]=lme(f,random=r,data=data,control=control,method=method,...)
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
	return(results)
}

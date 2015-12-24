# uses zcweights.environ and creates zcweights.environ.abun with environmental and fish abundance data
#
# Fishy abundance + environment
#
if(!exists("zcweights.environ"))
	source(file.path(system.file(package="CIPinnipedAnalysis"),"ZC_Weight_Environment_Model.r"))
# select data with fish abundance
zcweights.environ.abun=zcweights.environ[!is.na(zcweights.environ$SardineBiomass),]
zcweights.environ.abun=droplevels(zcweights.environ.abun)
# set up models
fixed.f.sa=c(fixed.f,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ SardineBiomass"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ AnchovyBiomass"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ hakeAge01"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ AnchovyBiomass + SardineBiomass"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ SardineAnchovyBiomass"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ log(SardineBiomass)"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ log(AnchovyBiomass)"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ log(SardineAnchovyBiomass)"))))
# fit models
random.f=list(res.environ$best.r)
res.environ.abun=fitmixed(fixed.f.sa,random.f,data=zcweights.environ.abun) 
res.environ.abun=compute_AICc(res.environ.abun, length(table(zcweights.environ.abun$Year))*2,5,fixed.f.sa)


res.environ.abun$model.table[1,]

# fit best model with reml
zc.weight.model.environ.abun=lme(res.environ.abun$best.f,random=res.environ.abun$best.r,data=zcweights.environ.abun,control=lmeControl(opt="optim"),method="REML")
summary(zc.weight.model.environ.abun)
zc.weight.model=zc.weight.model.environ.abun

bootstrap.se=function(x,nreps,days=0)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		cat("Bootstrap ",i,"\n")
		xsamp=lapply(split(x,list(x$batch,x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		mod=try(lme(fixed=res.environ.abun$best.f,random=res.environ.abun$best.r,data=as.data.frame(xsamp),method="REML",control=lmeControl(opt="optim")))
		if(class(mod)!="try-error")
		{
			i=i+1
			xsamp$days=days
			pp=predict(mod,newdata=xsamp)
			pmat[i,]=as.vector(tapply(as.vector(pp),list(x$cohort,x$sex),mean))
		}
	}
	return(sqrt(apply(pmat,2,var)))
}

#################  1 Oct Predictions ####################
# use 100 reps to compute std error
stderrors=bootstrap.se(zcweights.environ.abun,nboot,days=0)


pp=zcweights.environ.abun
pp$days=0
pp1=predict(zc.weight.model,newdata=pp)

pp0=predict(zc.weight.model,newdata=pp,level=0)

pp0=tapply(as.vector(pp0),list(zcweights.environ.abun$cohort,zcweights.environ.abun$sex),mean)
pp1=tapply(as.vector(pp1),list(zcweights.environ.abun$cohort,zcweights.environ.abun$sex),mean)

female.averages=data.frame(fit=pp1[,1],se=stderrors[1:length(pp1[,1])])
male.averages=data.frame(fit=pp1[,2],se=stderrors[(length(pp1[,1])+1):(2*length(pp1[,1]))])

expected.female.averages=data.frame(fit=pp0[,1],se=stderrors[1:length(pp0[,1])])
expected.male.averages=data.frame(fit=pp0[,2],se=stderrors[(length(pp0[,1])+1):(2*length(pp0[,1]))])

select=ZCWeight.df$Year%in%(as.numeric(rownames(pp1))+1975)

ZCWeight.df$female.environ.abun.mean.fall=NA
ZCWeight.df$female.environ.abun.mean.fall.se=NA
ZCWeight.df$male.environ.abun.mean.fall=NA
ZCWeight.df$male.environ.abun.mean.fall.se=NA
ZCWeight.df$male.environ.abun.mean.fall.fixed=NA
ZCWeight.df$female.environ.abun.mean.fall.fixed=NA

ZCWeight.df$female.environ.abun.mean.fall[select]=female.averages$fit
ZCWeight.df$female.environ.abun.mean.fall.se[select]=female.averages$se
ZCWeight.df$male.environ.abun.mean.fall[select]=male.averages$fit
ZCWeight.df$male.environ.abun.mean.fall.se[select]=male.averages$se

ZCWeight.df$female.environ.abun.mean.fall.fixed[select]=expected.female.averages$fit
ZCWeight.df$male.environ.abun.mean.fall.fixed[select]=expected.male.averages$fit

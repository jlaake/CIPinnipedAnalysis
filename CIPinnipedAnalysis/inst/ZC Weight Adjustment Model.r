require(CIPinnipedAnalysis)
if(!exists("fdir"))fdir=""
if(!exists("nboot"))nboot=100

# use "" to use databases in Calcur installed package directory; use NULL to use default Databases directory J:/Master  or specify directory
#fdir=NULL
#################################################################################
# Cross-sectional analysis
#################################################################################
# get zc weight values from database
zcweights=get.zc.weights(fdir=fdir)
#
#  exclude brand eval weights and captures in April and only use SMI
#
zcweights.all=zcweights[zcweights$days>-30&zcweights$days<150&zcweights$sitecode=="SMI",]
#
# compute observed fall averages - using data within 1 Sept to 15 Nov
#
zcweights=zcweights.all[zcweights.all$days>=-31&zcweights.all$days<=45,]
female.observed=sapply(split(zcweights$weight[zcweights$sex=="F"],
				factor(zcweights$cohort[zcweights$sex=="F"],levels=levels(factor(zcweights.all$cohort)))),mean)
male.observed=sapply(split(zcweights$weight[zcweights$sex=="M"],
				factor(zcweights$cohort[zcweights$sex=="M"],levels=levels(factor(zcweights.all$cohort)))),mean)

fall.avg.day=tapply(zcweights$days,zcweights$cohort,mean)
all.cohorts=sort(unique(zcweights.all$cohort))
#
# compute observed winter averages - using data > 15 Nov
#
zcweights=zcweights.all[zcweights.all$days>45,]
female.winter.observed=sapply(split(zcweights$weight[zcweights$sex=="F"],
				factor(zcweights$cohort[zcweights$sex=="F"],levels=levels(factor(zcweights.all$cohort)))),mean)
male.winter.observed=sapply(split(zcweights$weight[zcweights$sex=="M"],
				factor(zcweights$cohort[zcweights$sex=="M"],levels=levels(factor(zcweights.all$cohort)))),mean)
male.winter.observed[is.nan(male.winter.observed)]=NA
female.winter.observed[is.nan(female.winter.observed)]=NA

winter.avg.day=tapply(zcweights$days,factor(zcweights$cohort,levels=sort(unique(zcweights.all$cohort))),mean)
# For modelling use all of the data from 1 Sept to 1 Mar
zcweights=zcweights.all
zcweights$batch=factor(paste(zcweights$cohort,zcweights$days))
#
#  fit growth model- only adjusting for timing of collection date in fall
#
# First fit a sequence of random effect models with REML to assess best random model with same fixed model
random.f=list(list(~1|cohort,~1|batch),list(~days|cohort,~1|batch),list(~sex:days|cohort,~1|batch),list(~sex*days|cohort,~1|batch),list(~-1+sex+days|cohort,~1|batch))
fixed.f=list(weight~sex*days)
res.adjust=fitmixed(fixed.f,random.f,data=zcweights) 

# Using that random model, fit a sequence of fixed effect models with ML and use AIC to assess best fixed model
random.f=list(res.adjust$best.r)
fixed.f=list(
		weight~sex*days,
		weight~sex+days,
		weight~sex,
		weight~days,
		weight~1)
res.adjust=fitmixed(fixed.f,random.f,data=zcweights) 

# Finally fit best fixed/random model with REML
zc.weight.model=lme(fixed=res.adjust$best.f,random=res.adjust$best.r,data=zcweights,method="REML",control=lmeControl(opt="optim"))
print(summary(zc.weight.model))

# define bootstrap function to compute std error for predicted sex-cohort means
bootstrap.se=function(x,nreps,days)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		cat("Bootstrap ",i,"\n")
		xsamp=lapply(split(x,list(x$batch,x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		mod=try(lme(res.adjust$best.f,random=res.adjust$best.r,data=as.data.frame(xsamp),control=lmeControl(opt="optim")))
		if(class(mod)!="try-error")
		{
			i=i+1
			xsamp=data.frame(days=days,cohort=rep(sort(unique(x$cohort)),2),sex=rep(c("F","M"),each=length(unique(x$cohort))))
			pp=predict(mod,newdata=xsamp,level=1)
			pmat[i,]=unique(data.frame(predict=as.vector(pp),cohort=xsamp$cohort,sex=xsamp$sex))$predict
		}
	}
	return(sqrt(apply(pmat,2,var)))
}


################# 1 Oct Predictions #####################
# use nboot reps to compute std error
stderrors=bootstrap.se(zcweights,nboot,days=0)
# Compute predictions and construct dataframes for female and male averages with bootstrap std errors
pp=data.frame(days=0,cohort=rep(sort(unique(zcweights$cohort)),2),sex=rep(c("F","M"),each=length(unique(zcweights$cohort))))
pp$predict=predict(zc.weight.model,newdata=pp,level=1)
female.averages=data.frame(fit=pp$predict[pp$sex=="F"],se=stderrors[as.numeric(row.names(pp[pp$sex=="F",]))])
male.averages=data.frame(fit=pp$predict[pp$sex=="M"],se=stderrors[as.numeric(row.names(pp[pp$sex=="M",]))])
# create dataframe 
ZCWeight.df=data.frame(female.observed.mean.fall=female.observed,
		female.adjusted.mean.fall=female.averages$fit,female.adjusted.mean.fall.se=female.averages$se,
		male.observed.mean.fall=male.observed,
		male.adjusted.mean.fall=male.averages$fit,male.adjusted.mean.fall.se=male.averages$se,fall.avg.day=fall.avg.day)
ZCWeight.df=cbind(Year=as.numeric(rownames(ZCWeight.df)),ZCWeight.df)

options(width=150)
print(ZCWeight.df)
#
# Plot predicted female and male averages
#
pdf("ZCPredictedWeights.pdf",pointsize=10)
par(mfrow=c(2,1))
maxyear=max(zcweights$cohort)
ymin=min(c(female.averages$fit-1.96*female.averages$se, male.averages$fit-1.96*male.averages$se))*.9
ymax=max(c(female.averages$fit+1.96*female.averages$se, male.averages$fit+1.96*male.averages$se))*1.1
plot.weight.series(1975:maxyear,female.averages,main="California sea lion Female Pups",ylim=c(ymin,ymax),date="1 Oct")
plot.weight.series(1975:maxyear,male.averages,main="California sea lion Male Pups",ylim=c(ymin,ymax),date="1 Oct")
dev.off()

jpeg("ZCPredictedWeightsCalcofi.jpg",,height=600,width=600,quality=100,pointsize=12)
maxyear=max(zcweights$cohort)
par(lty=1)
plot.weight.series(1997:maxyear,female.averages[23:(maxyear-1974),],ylim=c(ymin,ymax),xaxp=c(1998,2012,7),date="1 Oct")
abline(h=mean(female.averages$fit[23:(maxyear-1974)]))
par(lty=2)
plot.weight.series(1997:maxyear,male.averages[23:(maxyear-1974),],pch=2,add=TRUE,slty=1,date="1 Oct")
abline(h=mean(male.averages$fit[23:(maxyear-1974)]))
points(2009,25,pch=2)
lines(x=c(2008.75,2009.25),y=c(25,25),pch=2,lty=2)
points(2009,23,pch=1)
lines(x=c(2008.75,2009.25),y=c(23,23),pch=1,lty=1)
text(2009.4,25,"Males",pos=4)
text(2009.4,23,"Females",pos=4)
par(lty=1)
dev.off()

################# 1 Feb Predictions #####################
pp=data.frame(days=123,cohort=rep(sort(unique(zcweights$cohort)),2),sex=rep(c("F","M"),each=length(unique(zcweights$cohort))))
pp$predict=predict(zc.weight.model,newdata=pp,level=1)
stderrors=bootstrap.se(zcweights,nboot,days=123)
female.averages=data.frame(fit=pp$predict[pp$sex=="F"],se=stderrors[as.numeric(row.names(pp[pp$sex=="F",]))])
male.averages=data.frame(fit=pp$predict[pp$sex=="M"],se=stderrors[as.numeric(row.names(pp[pp$sex=="M",]))])
ymin=min(c(female.averages$fit-1.96*female.averages$se, male.averages$fit-1.96*male.averages$se))*.9
ymax=max(c(female.averages$fit+1.96*female.averages$se, male.averages$fit+1.96*male.averages$se))*1.1

# add to dataframe
ZCWeight.df$female.observed.mean.winter=female.winter.observed
ZCWeight.df$female.adjusted.mean.winter=female.averages$fit
ZCWeight.df$female.adjusted.mean.se.winter=female.averages$se
ZCWeight.df$male.observed.mean.winter=male.winter.observed
ZCWeight.df$male.adjusted.mean.winter=male.averages$fit
ZCWeight.df$male.adjusted.mean.se.winter=male.averages$se
ZCWeight.df$winter.avg.day=winter.avg.day

jpeg("ZCPredictedWeightsFebCalcofi.jpg",height=600,width=600,quality=100,pointsize=12)
maxyear=max(zcweights$cohort)
par(lty=1)
plot.weight.series(1997:maxyear,female.averages[23:(maxyear-1974),],ylim=c(ymin,ymax),xaxp=c(1998,2012,7),date="1 Feb")
abline(h=mean(female.averages$fit[23:(maxyear-1974)]))
par(lty=2)
plot.weight.series(1997:maxyear,male.averages[23:(maxyear-1974),],pch=2,add=TRUE,slty=1,date="1 Feb")
abline(h=mean(male.averages$fit[23:(maxyear-1974)]))
points(2009,40,pch=2)
lines(x=c(2008.75,2009.25),y=c(40,40),pch=2,lty=2)
points(2009,38,pch=1)
lines(x=c(2008.75,2009.25),y=c(38,38),pch=1,lty=1)
text(2009.4,40,"Males",pos=4)
text(2009.4,38,"Females",pos=4)
par(lty=1)
dev.off()

# show plot of observed and predicted from adjustment mean for fall
jpeg("ZCAdjustment-ObservedComparison.jpg",height=600,width=600,quality=100,pointsize=12)
par(mfrow=c(2,1))
with(ZCWeight.df,plot(female.observed.mean.fall[fall.avg.day< -7],female.adjusted.mean.fall[fall.avg.day< -7],col="blue",ylim=c(12,23),xlim=c(12,23),xlab="Fall observed mean",ylab="Fall adjusted mean"))
with(ZCWeight.df,points(female.observed.mean.fall[fall.avg.day >= 7],female.adjusted.mean.fall[fall.avg.day>= 7],col="red"))
with(ZCWeight.df,points(female.observed.mean.fall[fall.avg.day >= -7 &fall.avg.day< 7],female.adjusted.mean.fall[fall.avg.day>= -7 &fall.avg.day<7],col="black"))
abline(0,1)
# show plot of observed and predicted from adjustment mean for fall
with(ZCWeight.df,plot(female.observed.mean.winter[winter.avg.day< 116],female.adjusted.mean.winter[winter.avg.day< 116],col="blue",ylim=c(14,32),xlim=c(14,32),xlab="Winter observed mean",ylab="Winter adjusted mean"))
with(ZCWeight.df,points(female.observed.mean.winter[winter.avg.day >= 130],female.adjusted.mean.winter[winter.avg.day>= 130],col="red"))
with(ZCWeight.df,points(female.observed.mean.winter[winter.avg.day >= 116 &winter.avg.day< 130],female.adjusted.mean.winter[winter.avg.day>= 116 &winter.avg.day< 130],col="black"))
abline(0,1)
dev.off()
# write out file for export to CIPinnipedCensusQuery
write.csv(ZCWeight.df,file="Zc weights.csv")


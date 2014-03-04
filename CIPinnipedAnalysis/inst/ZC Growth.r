require(CIPinnipedAnalysis)
if(!exists("fdir"))fdir=""
if(!exists("nboot"))nboot=100
# use "" to use databases in Calcur installed package directory; use NULL to use default Databases directory J:/Master or specify directory
#fdir=NULL
if(!exists("anomalies"))
{
	sdir=system.file(package="CIPinnipedAnalysis")
	source(file.path(sdir,"CreateAnomalies.r"))
}
################################################################################
# Longitudinal Environmental Growth Analysis
################################################################################
# get zc weight values from database
zcweights.lon=get.zc.weights(fdir=fdir)
#
#  Use SMI from 1 Sept to end of Feb
#
zcweights.lon=zcweights.lon[zcweights.lon$days>-30&zcweights.lon$days<150&zcweights.lon$sitecode=="SMI",]
# Extract the pups with more than one measurement
long.id=table(zcweights.lon$AnimalID)
long.id=names(long.id[long.id>1])
zcweights.lon=zcweights.lon[zcweights.lon$AnimalID%in%long.id,]
# Create initial and recap dataframes and exclude any with more than one initial
initial=zcweights.lon[zcweights.lon$type=="Initial",]
recaps=zcweights.lon[zcweights.lon$type=="Recap",]
long.id=table(initial$AnimalID)
long.id=names(long.id[long.id==1])
initial=initial[initial$AnimalID%in%long.id,]
zcweights.lon=rbind(initial,recaps)
# Compute growth rates
gr=sapply(split(zcweights.lon,list(zcweights.lon$AnimalID)),function(x) (diff(x$weight[order(x$days)])/diff(x$days[order(x$days)]))[1])

# Merge with data about animal
grdata=data.frame(gr=gr,AnimalID=names(gr))
grdata=merge(initial,grdata)

gr.mod=lme(gr~sex,random=list(~sex|cohort,~1|AnimalID),data=grdata)
pp=data.frame(cohort=rep(sort(unique(zcweights.lon$cohort)),2),sex=rep(c("F","M"),each=length(unique(zcweights.lon$cohort))))
pp$gr=predict(gr.mod,pp,level=1)

bootstrap.se=function(x,nreps)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		cat("Bootstrap ",i,"\n")
		xsamp=lapply(split(x,list(x$cohort,x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		xsamp$AnimalID=1:nrow(xsamp)
		mod=try(lme(formula(gr.mod),random=as.formula(gr.mod$call$random),data=as.data.frame(xsamp),control=lmeControl(opt="optim")))
		if(class(mod)!="try-error")
		{
			i=i+1
			pp=predict(mod,newdata=xsamp,level=1)
			pmat[i,]=as.vector(tapply(as.vector(pp),list(x$cohort,x$sex),mean))
		}
	}
	return(sqrt(apply(pmat,2,var)))
}
stderrors=bootstrap.se(grdata,nboot)
female.averages=data.frame(fit=pp$gr[pp$sex=="F"&pp$cohort>=1997],se=stderrors[as.numeric(row.names(pp[pp$sex=="F"&pp$cohort>=1997,]))])
male.averages=data.frame(fit=pp$gr[pp$sex=="M"&pp$cohort>=1997],se=stderrors[as.numeric(row.names(pp[pp$sex=="M"&pp$cohort>=1997,]))])

plot.growth.series=function (time, predictions,...)
{
	plotCI(time, predictions$fit, 1.96 * predictions$se,
			1.96 * predictions$se, xlab = "Cohort",
			ylab = paste("Average daily growth rate(kg/day)"),
			...)
	lines(time, predictions$fit)
	invisible()
}

jpeg("growth.jpg")
par(lty=1)
plot.growth.series(sort(unique(pp$cohort[pp$cohort>=1997])),female.averages,ylim=c(0,.12),xaxp=c(1998,2014,8))
par(lty=2)
plot.weight.series(sort(unique(pp$cohort[pp$cohort>=1997]))+.2,male.averages,pch=2,add=TRUE,slty=1,date="1 Oct")
points(2005,0.04,pch=2)
lines(x=c(2004.75,2005.25),y=c(.040,.040),pch=2,lty=2)
points(2005,0.03,pch=1)
lines(x=c(2004.75,2005.25),y=c(0.03,0.03),pch=1,lty=1)
text(2005.4,0.04,"Males",pos=4)
text(2005.4,0.03,"Females",pos=4)
dev.off()


# Merge growth data with environmental data
zcweights.environ.lon=merge(grdata,data.frame(cohort=1975:lastyear,SST=OcttoFebAnomalies[4:numyears],MEI=LaggedMEIOcttoFeb[-1],
				UWI33=UWImeansOcttoFeb[1,-(1:6)],UWI36=UWImeansOcttoFeb[2,-(1:6)]))
# Create model

fixed.f=list(gr~sex*UWI33,gr~sex*MEI,gr~sex*SST)
random.f=list(~1|cohort)

gr.models=fitmixed(fixed.f=fixed.f,random.f=random.f,data=zcweights.environ.lon)

mod=lme(fixed=gr.models$best.f,random=random.f,data=zcweights.environ.lon,control=)

plot(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="F"],predict(mod)[zcweights.environ.lon$sex=="F"],pch="F",type="b",ylim=c(0,.12))
points(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],predict(mod)[zcweights.environ.lon$sex=="M"],pch="M",type="b")
lines(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],predict(mod)[zcweights.environ.lon$sex=="M"],lty=2)
# Observed and predicted values
win.graph()
par(mfrow=c(1,2))
plot(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="F"],predict(mod)[zcweights.environ.lon$sex=="F"],pch="F",type="b",ylim=c(0,.14),ylab="Female daily growth rate (kg/day)",xlab="Cohort")
points(c(1989,1990,1995,1997:2010,2012:max(zcweights.environ.lon$cohort)),tapply(zcweights.environ.lon$gr[zcweights.environ.lon$sex=="F"],zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="F"],mean),pch="O")
plot(zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],predict(mod)[zcweights.environ.lon$sex=="M"],pch="M",type="b",ylim=c(0,.14),ylab="Male daily growth rate (kg/day)",xlab="Cohort")
points(c(1989,1990,1995,1997:2010,2012:max(zcweights.environ.lon$cohort)),tapply(zcweights.environ.lon$gr[zcweights.environ.lon$sex=="M"],zcweights.environ.lon$cohort[zcweights.environ.lon$sex=="M"],mean),pch="O")

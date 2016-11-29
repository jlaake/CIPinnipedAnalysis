#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
require(CIPinnipedAnalysis)
if(!exists("fdir"))fdir=NULL
if(!exists("nboot"))nboot=100
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
if(!exists("lastyear"))lastyear=max(zcweights.lon$cohort)
zcweights.lon=zcweights.lon[zcweights.lon$cohort<=lastyear,]

#
#  Use SMI from 1 Sept to end of Feb
#
zcweights.lon=zcweights.lon[zcweights.lon$days>-30&zcweights.lon$days<150&zcweights.lon$sitecode=="SMI",]
# to exclude Dec 2014 weights use
#zcweights.lon=zcweights.lon[zcweights.lon$days<70 | zcweights.lon$days>=78,]
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
growth= function(x){
	x=x[order(x$days),]
	(x$weight[length(x$weight)]-x$weight[1])/(x$days[length(x$days)]-x$days[1])	
} 
gr=sapply(split(zcweights.lon,list(zcweights.lon$AnimalID)),growth)

# Merge with data about animal
grdata=data.frame(gr=gr,AnimalID=names(gr))
grdata=merge(initial,grdata)

# to avoid duplicate records remove records with gr=Nan because days=0
grdata=grdata[!is.nan(grdata$gr),]
gr.mod=lme(gr~sex,random=list(~sex|cohort,~1|AnimalID),data=grdata,control=lmeControl(opt="optim"))
pp=data.frame(cohort=rep(sort(unique(zcweights.lon$cohort)),2),sex=rep(c("F","M"),each=length(unique(zcweights.lon$cohort))))
pp$gr=predict(gr.mod,pp,level=1)
pp$n=as.vector(table(grdata$cohort,grdata$sex))
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
pp$gr.se=stderrors

# Merge growth data with environmental data
if(length(1975:lastyear)>length(UWImeansOcttoFeb[1,-(1:6)]) & length(1975:lastyear)>length(OcttoFebAnomalies[-(1:3)]) & length(1975:lastyear)>length(LaggedMEIOcttoFeb[-1]) ){
	stop("\nmismatch in environmental data; specified lastyear is too great\n ")
} else {
	numyears=lastyear-1975+1
	zcweights.environ.lon=merge(grdata,data.frame(cohort=1975:lastyear,SST=OcttoFebAnomalies[-(1:3)][1:numyears],MEI=LaggedMEIOcttoFeb[-1][1:numyears],
					UWI33=UWImeansOcttoFeb[1,-(1:6)][1:numyears],UWI36=UWImeansOcttoFeb[2,-(1:6)][1:numyears]))
# Create model 
	fixed.f=list(gr~sex*UWI33,gr~sex*MEI,gr~sex*SST)
	random.f=list(~1|cohort)
	zcweights.environ.lon=subset(zcweights.environ.lon,select=c("cohort","gr","sex","UWI33","MEI","SST"))
	zcweights.environ.lon=na.exclude(zcweights.environ.lon)
	gr.models=fitmixed(fixed.f=fixed.f,random.f=random.f,data=zcweights.environ.lon)
	
	mod=lme(fixed=gr.models$best.f,random=random.f,data=zcweights.environ.lon)
} 

ppgr=data.frame(cohort=sort(unique(zcweights.lon$cohort)),sex=rep("F",each=length(unique(zcweights.lon$cohort))))
ppgr$gr=predict(gr.mod,ppgr,level=1)

pp.lon=predict(mod,type="response",level=0)
SST=sapply(split(zcweights.environ.lon$SST,zcweights.environ.lon$cohort),unique)
pred.gr=sapply(split(pp.lon[zcweights.environ.lon$sex=="F"],names(pp.lon[zcweights.environ.lon$sex=="F"])),unique)



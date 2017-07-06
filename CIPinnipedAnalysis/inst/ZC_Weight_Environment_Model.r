#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
require(CIPinnipedAnalysis)
if(!exists("nboot"))nboot=100
# If ZC_Weight_Adjustment_model has not been run, run now
if(!exists("ZCWeight.df"))
	source(file.path(system.file(package="CIPinnipedAnalysis"),"ZC_Weight_Adjustment_Model.r"))
# If anomalies don't exist create them now to use here
if(!exists("anomalies"))
	source(file.path(system.file(package="CIPinnipedAnalysis"),"CreateAnomalies.r"))
if(!exists("use.calcofi"))use.calcofi=FALSE
#################################################################################
# Cross-sectional analysis
#################################################################################
# get zc weight values from database
zcweights=get.zc.weights(fdir=fdir)
if(!exists("lastyear"))lastyear=max(zcweights$cohort)
zcweights=zcweights[zcweights$cohort<=lastyear,]
#
#  exclude brand eval weights and captures in April and only use SMI
#
zcweights=zcweights[zcweights$days>-30&zcweights$days<150&
				zcweights$sitecode=="SMI",]
zcweights$batch=factor(paste(zcweights$cohort,zcweights$days))
#
# Get environment data, fish abundance data and diet data
#
# Get calcofi and std environmental data
data(calcofi)
julycalcofi=calcofi[calcofi$Month=="July",]
julycalcofi=sapply(julycalcofi[julycalcofi$Station%in%c("76.7.49","76.7.51","76.7.55","80.51","80.55","83.3.51"),-(1:3)],function(x) tapply(x,julycalcofi$Year[julycalcofi$Station%in%c("76.7.49","76.7.51","76.7.55","80.51","80.55","83.3.51")],mean))
july.calcofi=t(t(julycalcofi)-colMeans(julycalcofi))
Octcalcofi=calcofi[calcofi$Month=="October",]
Octcalcofi=sapply(Octcalcofi[Octcalcofi$Station%in%c("76.7.49","76.7.51","76.7.55","80.51","80.55","83.3.51"),-(1:3)],function(x) tapply(x,Octcalcofi$Year[Octcalcofi$Station%in%c("76.7.49","76.7.51","76.7.55","80.51","80.55","83.3.51")],mean))
oct.calcofi=t(t(Octcalcofi)-colMeans(Octcalcofi))
colnames(oct.calcofi)=paste(colnames(oct.calcofi),".oct",sep="")
calcofi.df=as.data.frame(cbind(july.calcofi,oct.calcofi))
calcofi.df$cohort=as.numeric(rownames(calcofi.df))
numyears=lastyear-1975+1
if(length(1975:lastyear)>length(OcttoFebAnomalies[-(1:3)])) stop("lastyear set to value that is too great")
env.data=data.frame(cohort=1975:lastyear,SST=AprtoSeptAnomalies[-(1:3)][1:numyears],SST1=OcttoFebAnomalies[-(1:3)][1:numyears],SST2=JunetoFebAnomalies[-(1:3)][1:numyears],
		MEI=LaggedMEIAprtoSept[-1][1:numyears],MEI1=LaggedMEIOcttoFeb[-1][1:numyears],MEI2=LaggedMEIJunetoFeb[-1][1:numyears],UWI33=UWImeansAprtoSept[1,-(1:6)][1:numyears],UWI36=UWImeansAprtoSept[2,-(1:6)][1:numyears],
		UWI331=UWImeansOcttoFeb[1,-(1:6)][1:numyears],UWI361=UWImeansOcttoFeb[2,-(1:6)][1:numyears],UWI332=UWImeansJunetoFeb[1,-(1:6)][1:numyears],UWI362=UWImeansJunetoFeb[2,-(1:6)][1:numyears],
		SLH=AprtoSeptSLH[1:numyears],SLH1=OcttoFebSLH[1:numyears])
env.data=merge(env.data,calcofi.df,all.x=TRUE)

data(DensityPDONPGO)
env.data=cbind(env.data,DensityPDONPGO)
# get fish data
# read in abundance data files
data(sardine)
data(hake)
data(anchovy)
# compute age 0 to 1 hake, 0 to 2 and 1 to 2 hake
hake$Age01=hake$Age0+hake$Age1
hake$Age02=hake$Age0+hake$Age1+hake$Age2
hake$Age12=hake$Age1+hake$Age2
names(hake)[-1]=paste("hake",names(hake)[-1],sep="")
# add upper bound estimate for 2012-2013 anchovy
#anchovy=rbind(anchovy,data.frame(Year=2012:2013,AnchovyBiomass=c(30,30)))
# create sardine age 0-1 biomass
total.sa=data.frame(Year=1981:2013,Age01=sardine$Abundance[sardine$Age==0]+sardine$Abundance[sardine$Age==1],Age0=sardine$Abundance[sardine$Age==0],Age1=sardine$Abundance[sardine$Age==1])
total.sa$SardineBiomass=sardine$Biomass[sardine$Age==1]+sardine$Biomass[sardine$Age==0]
names(total.sa)[2:4]=paste("sardine",names(total.sa)[2:4],sep="")
# convert abundance to 1000s and biomass to 1000mt
total.sa[,2:5]=total.sa[,2:5]/1000
fish_abundance=merge(total.sa,anchovy)
# compute biomass of anchovy and sardines in 1000s metric tons; excludes some years when only one of the species 
fish_abundance$SardineAnchovyBiomass=fish_abundance$SardineBiomass+fish_abundance$AnchovyBiomass
#merge hake data
fish_abundance=merge(hake,fish_abundance,all.x=TRUE)
names(fish_abundance)[1]="cohort"
#merge fish abundance and environmental data
env.data=merge(env.data,fish_abundance,all.x=TRUE )
# get diet data
# create freq of occurence for prey data from dataframe fo in
# package
data(fo)
fo=create_fo(fo)
fo$diet=NA
fo$diet[fo$Year%in%1981:1986]="Diet1"
fo$diet[fo$Year%in%c(1980,2009,2012,2013,1992,2000,1995,2001,2011,1991,2010,2014,2015)]="Diet2"
fo$diet[fo$Year%in%c(1993,1996,1997,1998,2002,2003,2004,2005,2006,2007)]="Diet3"
fo$diet=factor(fo$diet)
names(fo)[1]="cohort"
#merge in diet data
env.data=merge(env.data,fo,all.x=TRUE)
names(env.data)[1]="Year"
if(any(sapply(env.data,function(x) any(is.na(x)))))
{
   cat("Following variables are missing data\n")
   cat(paste(names(env.data)[sapply(env.data,function(x) any(is.na(x)))],sep="\n"))
}
#
# merge weight data with environment data
zcweights.environ=merge(zcweights,env.data,all.x=TRUE,by.x="cohort",by.y="Year")
zcweights.environ$Year=zcweights.environ$cohort
zcweights.environ$cohort=zcweights.environ$Year-min(zcweights.environ$Year)
# limit data to calcofi years if use.calcofi set to TRUE
if(use.calcofi)
{
	zcweights.environ=zcweights.environ[!is.na(zcweights.environ$dynamic_height_0_500m),]
	zcweights.environ=droplevels(zcweights.environ)
}	
#
# specify set of models
source(file.path(system.file(package="CIPinnipedAnalysis"),"environment_models.r"))
#
# First fit a sequence of random effect models with REML to assess best random model with same fixed model
random.f=list(list(~1|cohort,~1|batch),list(~days|cohort,~1|batch),list(~sex:days|cohort,~1|batch),list(~sex*days|cohort,~1|batch),list(~-1+sex+days|cohort,~1|batch))
fixed.f1=list(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+UWI33+cohort+MEI)
res.environ=fitmixed(fixed.f1,random.f,data=zcweights.environ) 

# Using that random model, fit a sequence of fixed effect models with ML and use AIC to assess best fixed model
random.f=list(res.environ$best.r)
res.environ=fitmixed(fixed.f,random.f,data=zcweights.environ) 
res.environ=compute_AICc(res.environ, length(table(zcweights.environ$Year))*2,5,fixed.f)

# Finally fit best fixed/random model with REML
zc.weight.environ.model=lme(fixed=res.environ$best.f,random=res.environ$best.r,data=res.environ$data,method="REML",control=lmeControl(opt="optim"))
print(summary(zc.weight.environ.model))


bootstrap.se=function(x,nreps,days=0)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		cat("Bootstrap ",i,"\n")
		xsamp=lapply(split(x,list(x$batch,x$sex)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		mod=try(lme(fixed=res.environ$best.f,random=res.environ$best.r,data=as.data.frame(xsamp),method="REML",control=lmeControl(opt="optim")))
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
stderrors=bootstrap.se(res.environ$data,nboot,days=0)
# Compute fall 1 Oct predictions and construct dataframes for female and male averages with std errors
pp=res.environ$data
pp$days=0
# predictions at 1 Oct with random effects
pp1=predict(zc.weight.environ.model,newdata=pp)
# predictions at 1 Oct with fixed effects only
pp0=predict(zc.weight.environ.model,newdata=pp,level=0)
# compute mean values which essentially acts as unique
pp0=tapply(as.vector(pp0),list(pp$cohort,pp$sex),mean)
pp1=tapply(as.vector(pp1),list(pp$cohort,pp$sex),mean)
# create list with estimates and std errors for averages
female.averages=data.frame(fit=pp1[,1],se=stderrors[1:length(pp1[,1])])
male.averages=data.frame(fit=pp1[,2],se=stderrors[(length(pp1[,1])+1):(2*length(pp1[,1]))])
# create list with estimates and std errors for expected averages
expected.female.averages=data.frame(fit=pp0[,1],se=stderrors[1:length(pp0[,1])])
expected.male.averages=data.frame(fit=pp0[,2],se=stderrors[(length(pp0[,1])+1):(2*length(pp0[,1]))])

ZCWeight.df$female.environ.mean.fall=NA
ZCWeight.df$female.environ.mean.fall.se=NA
ZCWeight.df$male.environ.mean.fall=NA
ZCWeight.df$male.environ.mean.fall.se=NA
ZCWeight.df$female.environ.mean.fall.fixed=NA
ZCWeight.df$male.environ.mean.fall.fixed=NA

add=0
if(use.calcofi)add=9
ZCWeight.df$female.environ.mean.fall[(add+1):(add+length(female.averages$fit))]=female.averages$fit
ZCWeight.df$female.environ.mean.fall.se[(add+1):(add+length(female.averages$fit))]=female.averages$se
ZCWeight.df$male.environ.mean.fall[(add+1):(add+length(female.averages$fit))]=male.averages$fit
ZCWeight.df$male.environ.mean.fall.se[(add+1):(add+length(female.averages$fit))]=male.averages$se
	
ZCWeight.df$female.environ.mean.fall.fixed[(add+1):(add+length(female.averages$fit))]=expected.female.averages$fit
ZCWeight.df$male.environ.mean.fall.fixed[(add+1):(add+length(female.averages$fit))]=expected.male.averages$fit
	

#################  1 Feb Predictions ####################
# use 100 reps to compute std error
stderrors=bootstrap.se(res.environ$data,nboot,days=123)
# Compute fall 1 Feb predictions and construct dataframes for female and male averages with std errors
pp=res.environ$data
pp$days=123
# predictions at 1 Oct with random effects
pp1=predict(zc.weight.environ.model,newdata=pp)
# predictions at 1 Oct with fixed effects only
pp0=predict(zc.weight.environ.model,newdata=pp,level=0)
# compute mean values which essentially acts as unique
pp0=tapply(as.vector(pp0),list(pp$cohort,pp$sex),mean)
pp1=tapply(as.vector(pp1),list(pp$cohort,pp$sex),mean)
# create list with estimates and std errors for averages
female.averages=data.frame(fit=pp1[,1],se=stderrors[1:length(pp1[,1])])
male.averages=data.frame(fit=pp1[,2],se=stderrors[(length(pp1[,1])+1):(2*length(pp1[,1]))])
# create list with estimates and std errors for expected averages
expected.female.averages=data.frame(fit=pp0[,1],se=stderrors[1:length(pp0[,1])])
expected.male.averages=data.frame(fit=pp0[,2],se=stderrors[(length(pp0[,1])+1):(2*length(pp0[,1]))])

ZCWeight.df$female.environ.mean.winter=NA
ZCWeight.df$female.environ.mean.winter.se=NA
ZCWeight.df$female.environ.mean.winter.fixed=NA

ZCWeight.df$male.environ.mean.winter=NA
ZCWeight.df$male.environ.mean.winter.se=NA
ZCWeight.df$male.environ.mean.winter.fixed=NA
	
ZCWeight.df$female.environ.mean.winter[(add+1):(add+length(female.averages$fit))]=female.averages$fit
ZCWeight.df$female.environ.mean.winter.se[(add+1):(add+length(female.averages$fit))]=female.averages$se
ZCWeight.df$male.environ.mean.winter[(add+1):(add+length(female.averages$fit))]=male.averages$fit
ZCWeight.df$male.environ.mean.winter.se[add+1:length(female.averages$fit)]=male.averages$se

ZCWeight.df$female.environ.mean.winter.fixed[(add+1):(add+length(female.averages$fit))]=expected.female.averages$fit
ZCWeight.df$male.environ.mean.winter.fixed[(add+1):(add+length(female.averages$fit))]=expected.male.averages$fit


	
		
		
		
		

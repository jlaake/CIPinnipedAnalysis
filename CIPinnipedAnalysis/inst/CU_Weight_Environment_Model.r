#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
if(!exists("nboot"))nboot=100
if(!exists("use.calcofi"))use.calcofi=FALSE
if(!exists("CUWeight.df"))
	source(file.path(system.file(package="CIPinnipedAnalysis"),"CU_Weight_Adjustment_Model.r"))


################################################################################
# Construct environmental values for models
################################################################################
if(!exists("anomalies"))
{
	sdir=system.file(package="CIPinnipedAnalysis")
	source(file.path(sdir,"CreateAnomalies.r"))
}


# get cu weight values from database
cuweights=get.cu.weights(fdir=fdir)
if(!exists("lastyear"))lastyear=max(cuweights$cohort)
#
#  exclude brand eval weights and captures in April
#
cuweights=cuweights[cuweights$days>-30&cuweights$days<150,]
cuweights=cuweights[cuweights$cohort<=lastyear,]
#
# Get environment data, fish abundance data and diet data
#
# Get calcofi and std environmental data
data(calcofi)
calcofi=calcofi[calcofi$Month=="July",]
calcofi=sapply(calcofi[calcofi$Station%in%c("76.7.49","76.7.51","76.7.55","80.51","80.55","83.3.51"),-(1:3)],function(x) tapply(x,calcofi$Year[calcofi$Station%in%c("76.7.49","76.7.51","76.7.55","80.51","80.55","83.3.51")],mean))
july.calcofi=t(t(calcofi)-colMeans(calcofi))
data(calcofi)
calcofi=calcofi[calcofi$Month=="October",]
calcofi=sapply(calcofi[calcofi$Station%in%c("76.7.49","76.7.51","76.7.55","80.51","80.55","83.3.51"),-(1:3)],function(x) tapply(x,calcofi$Year[calcofi$Station%in%c("76.7.49","76.7.51","76.7.55","80.51","80.55","83.3.51")],mean))
oct.calcofi=t(t(calcofi)-colMeans(calcofi))
colnames(oct.calcofi)=paste(colnames(oct.calcofi),".oct",sep="")
calcofi=as.data.frame(cbind(july.calcofi,oct.calcofi))
calcofi$cohort=as.numeric(rownames(calcofi))
numyears=lastyear-1975+1
if(length(1975:lastyear)>length(OcttoFebAnomalies[-(1:3)])) stop("lastyear set to value that is too great")
env.data=data.frame(cohort=1975:lastyear,SST=AprtoSeptAnomalies[-(1:3)][1:numyears],SST1=OcttoFebAnomalies[-(1:3)][1:numyears],SST2=JunetoFebAnomalies[-(1:3)][1:numyears],
		MEI=LaggedMEIAprtoSept[-1][1:numyears],MEI1=LaggedMEIOcttoFeb[-1][1:numyears],MEI2=LaggedMEIJunetoFeb[-1][1:numyears],UWI33=UWImeansAprtoSept[1,-(1:6)][1:numyears],UWI36=UWImeansAprtoSept[2,-(1:6)][1:numyears],
		UWI331=UWImeansOcttoFeb[1,-(1:6)][1:numyears],UWI361=UWImeansOcttoFeb[2,-(1:6)][1:numyears],UWI332=UWImeansJunetoFeb[1,-(1:6)][1:numyears],UWI362=UWImeansJunetoFeb[2,-(1:6)][1:numyears],
		SLH=AprtoSeptSLH[1:numyears],SLH1=OcttoFebSLH[1:numyears])

env.data=merge(env.data,calcofi,all.x=TRUE)
# get fish data
# read in abundance data files
#sardine=read.delim("sardine.txt",header=TRUE)
#anchovy=read.delim("AnchovyAbundance.txt",header=TRUE)
#hake=read.delim("hake.txt",header=TRUE)
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

cuweights.environ=merge(cuweights,env.data,all.x=TRUE,by.x="cohort",by.y="cohort")
cuweights.environ$Year=cuweights.environ$cohort
cuweights.environ$cohort=cuweights.environ$Year-min(cuweights.environ$Year)
# limit data to calcofi years if use.calcofi set to TRUE
if(use.calcofi)
{
	cuweights.environ=cuweights.environ[!is.na(cuweights.environ$dynamic_height_0_500m),]
	cuweights.environ=droplevels(cuweights.environ)
}	
#cuweights.environ=merge(cuweights,data.frame(cohort=1975:lastyear,SST=JunetoSeptAnomalies[4:numyears],SST1=OcttoFebAnomalies[4:numyears],SST2=JunetoFebAnomalies[4:numyears],
#				MEI=LaggedMEIJunetoSept[-1],MEI1=LaggedMEIOcttoFeb[-1],MEI2=LaggedMEIJunetoFeb[-1],UWI33=UWImeansJunetoSept[1,-(1:6)],UWI36=UWImeansJunetoSept[2,-(1:6)],
#				UWI331=UWImeansOcttoFeb[1,-(1:6)],UWI361=UWImeansOcttoFeb[2,-(1:6)],UWI332=UWImeansJunetoFeb[1,-(1:6)],UWI362=UWImeansJunetoFeb[2,-(1:6)]))
#cuweights.environ$cohort.factor=factor(ifelse(cuweights.environ$cohort<1990,0,1),labels=c("<=1989",">=1990"))
#cuweights.environ$SST=cuweights.environ$SST1
#cuweights.environ$SST[cuweights.environ$days>90]=cuweights.environ$SST3[cuweights.environ$days>90]
#cuweights.environ$MEI[cuweights.environ$days>90]=cuweights.environ$MEI1[cuweights.environ$days>90]
#cuweights.environ$UWI33[cuweights.environ$days>90]=cuweights.environ$UWI331[cuweights.environ$days>90]
#cuweights.environ$UWI36[cuweights.environ$days>90]=cuweights.environ$UWI361[cuweights.environ$days>90]

#############################################################################################################
# Environmental Cross-sectional analysis of weights
#
# First fit a sequence of random effect models with REML to assess best random model with same fixed model
random.f=list(list(~1|cohort),list(~days|cohort),list(~sex:days|cohort),list(~sex*days|cohort),list(~-1+sex+days|cohort))
fixed.f=list(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+UWI33+cohort+MEI)
res.environ=fitmixed(fixed.f,random.f,data=cuweights.environ) 

# Using that random model, fit a sequence of fixed effect models with ML and use AIC to assess best fixed model
random.f=list(res.environ$best.r)
# specify set of models
source(file.path(system.file(package="CIPinnipedAnalysis"),"environment_models.r"))

res.environ=fitmixed(fixed.f,random.f,data=cuweights.environ) 
res.environ=compute_AICc(res.environ, length(table(cuweights.environ$Year))*2,4,fixed.f)

# Finally fit best fixed/random model with REML
cu.weight.environ.model=lme(fixed=res.environ$best.f,random=res.environ$best.r,data=res.environ$data,method="REML",control=lmeControl(opt="optim"))
print(summary(cu.weight.environ.model))

bootstrap.se=function(x,nreps)
{
	pmat=matrix(0,nrow=nreps,ncol=nrow(unique(data.frame(cohort=x$cohort,sex=x$sex))))
	i=0
	while(i<nreps)
	{
		cat("Bootstrap ",i,"\n")
		xsamp=lapply(split(x,list(x$sex,x$cohort)),function(x) if(nrow(x)>0) x[sample(1:nrow(x),replace=TRUE),] else NULL)
		xsamp=do.call("rbind",xsamp)
		mod=try(lme(formula(cu.weight.model),random=as.formula(cu.weight.model$call$random),data=as.data.frame(xsamp),control=lmeControl(opt="optim")))
		if(class(mod)!="try-error")
		{
			i=i+1
			xsamp$days=0
			pp=predict(mod,newdata=xsamp)
			pmat[i,]=as.vector(tapply(as.vector(pp),list(x$cohort,x$sex),mean))
		}
	}
	return(sqrt(apply(pmat,2,var)))
}
# use 100 reps to compute std error
stderrors=bootstrap.se(res.environ$data,nboot)
# Compute fall 1 Oct predictions and construct dataframes for female and male averages with std errors
pp=res.environ$data
pp$days=0
# predictions at 1 Oct with random effects
pp1=predict(cu.weight.environ.model,newdata=pp)
# predictions at 1 Oct with fixed effects only
pp0=predict(cu.weight.environ.model,newdata=pp,level=0)
# compute mean values which essentially acts as unique
pp0=tapply(as.vector(pp0),list(pp$cohort,pp$sex),mean)
pp1=tapply(as.vector(pp1),list(pp$cohort,pp$sex),mean)
# create list with estimates and std errors for averages
female.averages=data.frame(fit=pp1[,1],se=stderrors[1:length(pp1[,1])])
male.averages=data.frame(fit=pp1[,2],se=stderrors[(length(pp1[,1])+1):(2*length(pp1[,1]))])
# create list with estimates and std errors for expected averages
expected.female.averages=data.frame(fit=pp0[,1],se=stderrors[1:length(pp0[,1])])
expected.male.averages=data.frame(fit=pp0[,2],se=stderrors[(length(pp0[,1])+1):(2*length(pp0[,1]))])

CUWeight.df$female.environ.mean.fall=NA
CUWeight.df$female.environ.mean.fall.se=NA
CUWeight.df$male.environ.mean.fall=NA
CUWeight.df$male.environ.mean.fall.se=NA
CUWeight.df$female.environ.mean.fall.fixed=NA
CUWeight.df$male.environ.mean.fall.fixed=NA

add=0
if(use.calcofi)add=9
CUWeight.df$female.environ.mean.fall[(add+1):(add+length(female.averages$fit))]=female.averages$fit
CUWeight.df$female.environ.mean.fall.se[(add+1):(add+length(female.averages$fit))]=female.averages$se
CUWeight.df$male.environ.mean.fall[(add+1):(add+length(female.averages$fit))]=male.averages$fit
CUWeight.df$male.environ.mean.fall.se[(add+1):(add+length(female.averages$fit))]=male.averages$se

CUWeight.df$female.environ.mean.fall.fixed[(add+1):(add+length(female.averages$fit))]=expected.female.averages$fit
CUWeight.df$male.environ.mean.fall.fixed[(add+1):(add+length(female.averages$fit))]=expected.male.averages$fit

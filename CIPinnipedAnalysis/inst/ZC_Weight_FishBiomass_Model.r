# uses zcweights.environ and creates zcweights.environ.abun with environmental and fish abundance data
#
# Sardine and anchovy abundance + environment
#
if(!exists("zcweights.environ"))
	source(file.path(system.file(package="CIPinnipedAnalysis"),"ZC_Weight_Environment_Model.r"))

# read in abundance data files
	sardine=read.delim("sardine.txt",header=T)
anchovy=read.delim("AnchovyAbundance.txt",header=T)
# add rought upper bound estimate for 2012-2013
anchovy=rbind(anchovy,data.frame(Year=2012:2013,AnchovyBiomass=c(20,20)))
total.sa=data.frame(Year=1981:2013,Age01=sardine$Abundance[sardine$Age==0]+sardine$Abundance[sardine$Age==1],Age0=sardine$Abundance[sardine$Age==0],Age1=sardine$Abundance[sardine$Age==1])
total.sa$SardineBiomass=sardine$Biomass[sardine$Age==1]+sardine$Biomass[sardine$Age==0]
# convert abundance to 1000s and biomass to 1000mt
total.sa[,2:5]=total.sa[,2:5]/1000
fish_abundance=merge(total.sa,anchovy)
# compute biomass of anchovy and sardines in 1000s metric tons; excludes some years when only one of the species
fish_abundance$SardineAnchovyBiomass=fish_abundance$SardineBiomass+fish_abundance$AnchovyBiomass
#
if(use.calcofi)
{
	zcweights.environ$Year=zcweights.environ$cohort+1984
	fixed.f=list(		
			weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort,
			weight~sex*SST+sex:days+SST1:days+sex:SST1:days,
			weight~sex*SST+sex:days+SST1:days+cohort,
			weight~sex*SST+sex:days+SST1:days,
			weight~sex*SST+sex:days+cohort,
			weight~sex*SST+sex:days,
			weight~sex+SST+sex:days+SST1:days+sex:SST1:days+cohort,
			weight~sex+SST+sex:days+SST1:days+sex:SST1:days,
			weight~sex+SST+sex:days+SST1:days+cohort,
			weight~sex+SST+sex:days+SST1:days,
			weight~sex+SST+sex:days+cohort,
			weight~sex+SST+sex:days,
			
			weight~sex*dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days+sex:dynamic_height_0_500m.oct:days+cohort,
			weight~sex*dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days+sex:dynamic_height_0_500m.oct:days,
			weight~sex*dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days+cohort,
			weight~sex*dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days,
			weight~sex*dynamic_height_0_500m+sex:days+cohort,
			weight~sex*dynamic_height_0_500m+sex:days,
			weight~sex+dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days+sex:dynamic_height_0_500m.oct:days+cohort,
			weight~sex+dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days+sex:dynamic_height_0_500m.oct:days,
			weight~sex+dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days+cohort,
			weight~sex+dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days,
			weight~sex+dynamic_height_0_500m+sex:days+cohort,
			weight~sex+dynamic_height_0_500m+sex:days,
			
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+sex:R_SIGMA_75m.oct:days+cohort,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+sex:R_SIGMA_75m.oct:days,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days,
			weight~sex*R_SIGMA_75m+sex:days+cohort,
			weight~sex*R_SIGMA_75m+sex:days,
			weight~sex+R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+sex:R_SIGMA_75m.oct:days+cohort,
			weight~sex+R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+sex:R_SIGMA_75m.oct:days,
			weight~sex+R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort,
			weight~sex+R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days,
			weight~sex+R_SIGMA_75m+sex:days+cohort,
			weight~sex+R_SIGMA_75m+sex:days,
			
			weight~sex*R_POTEMP_75m+sex:days+R_POTEMP_75m.oct:days+sex:R_POTEMP_75m.oct:days+cohort,
			weight~sex*R_POTEMP_75m+sex:days+R_POTEMP_75m.oct:days+sex:R_POTEMP_75m.oct:days,
			weight~sex*R_POTEMP_75m+sex:days+R_POTEMP_75m.oct:days+cohort,
			weight~sex*R_POTEMP_75m+sex:days+R_POTEMP_75m.oct:days,
			weight~sex*R_POTEMP_75m+sex:days+cohort,
			weight~sex*R_POTEMP_75m+sex:days,
			weight~sex+R_POTEMP_75m+sex:days+R_POTEMP_75m.oct:days+sex:R_POTEMP_75m.oct:days+cohort,
			weight~sex+R_POTEMP_75m+sex:days+R_POTEMP_75m.oct:days+sex:R_POTEMP_75m.oct:days,
			weight~sex+R_POTEMP_75m+sex:days+R_POTEMP_75m.oct:days+cohort,
			weight~sex+R_POTEMP_75m+sex:days+R_POTEMP_75m.oct:days,
			weight~sex+R_POTEMP_75m+sex:days+cohort,
			weight~sex+R_POTEMP_75m+sex:days,
			
			weight~sex*stratification+sex:days+stratification.oct:days+sex:stratification.oct:days+cohort,
			weight~sex*stratification+sex:days+stratification.oct:days+sex:stratification.oct:days,
			weight~sex*stratification+sex:days+stratification.oct:days+cohort,
			weight~sex*stratification+sex:days+stratification.oct:days,
			weight~sex*stratification+sex:days+cohort,
			weight~sex*stratification+sex:days,
			weight~sex+stratification+sex:days+stratification.oct:days+sex:stratification.oct:days+cohort,
			weight~sex+stratification+sex:days+stratification.oct:days+sex:stratification.oct:days,
			weight~sex+stratification+sex:days+stratification.oct:days+cohort,
			weight~sex+stratification+sex:days+stratification.oct:days,
			weight~sex+stratification+sex:days+cohort,
			weight~sex+stratification+sex:days,
			
			weight~sex*pycnocline_depth+sex:days+pycnocline_depth.oct:days+sex:pycnocline_depth.oct:days+cohort,
			weight~sex*pycnocline_depth+sex:days+pycnocline_depth.oct:days+sex:pycnocline_depth.oct:days,
			weight~sex*pycnocline_depth+sex:days+pycnocline_depth.oct:days+cohort,
			weight~sex*pycnocline_depth+sex:days+pycnocline_depth.oct:days,
			weight~sex*pycnocline_depth+sex:days+cohort,
			weight~sex*pycnocline_depth+sex:days,
			weight~sex+pycnocline_depth+sex:days+pycnocline_depth.oct:days+sex:pycnocline_depth.oct:days+cohort,
			weight~sex+pycnocline_depth+sex:days+pycnocline_depth.oct:days+sex:pycnocline_depth.oct:days,
			weight~sex+pycnocline_depth+sex:days+pycnocline_depth.oct:days+cohort,
			weight~sex+pycnocline_depth+sex:days+pycnocline_depth.oct:days,
			weight~sex+pycnocline_depth+sex:days+cohort,
			weight~sex+pycnocline_depth+sex:days)
			
} else{
    fixed.f=list(
		weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort,
		weight~sex*SST+sex:days+SST1:days+sex:SST1:days,
		weight~sex*SST+sex:days+SST1:days+cohort,
		weight~sex*SST+sex:days+SST1:days,
		weight~sex*SST+sex:days+cohort,
		weight~sex*SST+sex:days,
		weight~sex+SST+sex:days+SST1:days+sex:SST1:days+cohort,
		weight~sex+SST+sex:days+SST1:days+sex:SST1:days,
		weight~sex+SST+sex:days+SST1:days+cohort,
		weight~sex+SST+sex:days+SST1:days,
		weight~sex+SST+sex:days+cohort,
		weight~sex+SST+sex:days)
	zcweights.environ$Year=zcweights.environ$cohort+1975
}
zcweights.environ.abun=merge(zcweights.environ,fish_abundance,by="Year")

fixed.f.sa=c(fixed.f,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ SardineBiomass"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ AnchovyBiomass"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ SardineAnchovyBiomass"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ log(SardineBiomass)"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ log(AnchovyBiomass)"))))
fixed.f.sa=c(fixed.f.sa,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ log(SardineAnchovyBiomass)"))))

res.environ.abun=fitmixed(fixed.f.sa,random.f,data=zcweights.environ.abun) 

res.environ.abun$model.table[1,]


zc.weight.model.environ.abun=lme(res.environ.abun$best.f,random=res.environ.abun$best.r,data=zcweights.environ.abun,control=lmeControl(opt="optim"),method="REML")
summary(zc.weight.model.environ.abun)
zc.weight.model=zc.weight.model.environ.abun

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

if(use.calcofi){
	select=ZCWeight.df$Year%in%(as.numeric(rownames(pp1))+1984)
}else{
	select=ZCWeight.df$Year%in%(as.numeric(rownames(pp1))+1975)
}
ZCWeight.df$female.eviron.abun.mean.fall=NA
ZCWeight.df$female.eviron.abun.mean.fall.se=NA
ZCWeight.df$male.eviron.abun.mean.fall=NA
ZCWeight.df$male.eviron.abun.mean.fall.se=NA
ZCWeight.df$male.eviron.abun.mean.fall.fixed=NA
ZCWeight.df$female.eviron.abun.mean.fall.fixed=NA

ZCWeight.df$female.eviron.abun.mean.fall[select]=female.averages$fit
ZCWeight.df$female.eviron.abun.mean.fall.se[select]=female.averages$se
ZCWeight.df$male.eviron.abun.mean.fall[select]=male.averages$fit
ZCWeight.df$male.eviron.abun.mean.fall.se[select]=male.averages$se

ZCWeight.df$female.environ.abun.mean.fall.fixed[select]=expected.female.averages$fit
ZCWeight.df$male.environ.mean.fall.fixed[select]=expected.male.averages$fit

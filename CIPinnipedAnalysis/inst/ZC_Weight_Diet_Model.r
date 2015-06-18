#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
require(CIPinnipedAnalysis)
if(!exists("lastyear"))lastyear=2013
if(!exists("nboot"))nboot=100
if(!exists("anomalies"))
if(!exists("use.calcofi"))use.calcofi=FALSE
# if ZC_Weight_Environment_Model has not been run, run it now
if(!exists("ZCWeight.df"))
{
	source(file.path(sdir,"ZC_Weight_Environment_Model.r"))
} else 
{
	if(is.null(ZCWeight.df$female.environ.mean.fall))
		source(file.path(sdir,"ZC_Weight_Environment_Model.r"))
}
# set of models depends on whether calcofi data are used
if(use.calcofi)
{
	data(calcofi)
	calcofi=calcofi[calcofi$Month=="July",]
	calcofi=sapply(calcofi[,-(1:3)],function(x) tapply(x,calcofi$Year,mean))
	july.calcofi=t(t(calcofi)-colMeans(calcofi))
	data(calcofi)
	calcofi=calcofi[calcofi$Month=="October",]
	calcofi=sapply(calcofi[,-(1:3)],function(x) tapply(x,calcofi$Year,mean))
	oct.calcofi=t(t(calcofi)-colMeans(calcofi))
	colnames(oct.calcofi)=paste(colnames(oct.calcofi),".oct",sep="")
	zcweights.environ.diet=merge(zcweights,cbind(july.calcofi,oct.calcofi,data.frame(cohort=1984:lastyear,SST=JunetoSeptAnomalies[13:numyears],
							SST1=OcttoFebAnomalies[13:numyears],SST2=JunetoFebAnomalies[13:numyears],
							MEI=LaggedMEIJunetoSept[-(1:10)],MEI1=LaggedMEIOcttoFeb[-(1:10)],MEI2=LaggedMEIJunetoFeb[-(1:10)],
							UWI33=UWImeansJunetoSept[1,-(1:15)],UWI36=UWImeansJunetoSept[2,-(1:15)],
							UWI331=UWImeansOcttoFeb[1,-(1:15)],UWI361=UWImeansOcttoFeb[2,-(1:15)],UWI332=UWImeansJunetoFeb[1,-(1:15)],
							UWI362=UWImeansJunetoFeb[2,-(1:15)])))
	fixed.f=list(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort.factor,
			weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort,
			weight~sex*SST+sex:days+SST1:days+sex:SST1:days,
			weight~sex*SST+sex:days+SST1:days+cohort.factor,
			weight~sex*SST+sex:days+SST1:days+cohort,
			weight~sex*SST+sex:days+SST1:days,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days+cohort.factor,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days+cohort,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days,
			weight~sex*UWI36+sex:days+UWI361:days+cohort.factor,
			weight~sex*UWI36+sex:days+UWI361:days+cohort,
			weight~sex*UWI36+sex:days+UWI361:days,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days+cohort.factor,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days+cohort,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days,
			weight~sex*MEI+sex:days+MEI1:days+cohort.factor,
			weight~sex*MEI+sex:days+MEI1:days+cohort,
			weight~sex*MEI+sex:days,
			weight~sex*dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days+sex:dynamic_height_0_500m.oct:days+cohort.factor,
			weight~sex*dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days+sex:dynamic_height_0_500m.oct:days+cohort,
			weight~sex*dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days+sex:dynamic_height_0_500m.oct:days,
			weight~sex*dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days+cohort.factor,
			weight~sex*dynamic_height_0_500m+sex:days+dynamic_height_0_500m.oct:days+cohort,
			weight~sex*dynamic_height_0_500m+sex:days,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+sex:R_SIGMA_75m.oct:days+cohort.factor,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+sex:R_SIGMA_75m.oct:days+cohort,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+sex:R_SIGMA_75m.oct:days,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort.factor,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort,
			weight~sex*R_SIGMA_75m+sex:days,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m.oct:days+sex:R_POTEMP_25m.oct:days+cohort.factor,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m.oct:days+sex:R_POTEMP_25m.oct:days+cohort,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m.oct:days+sex:R_POTEMP_25m.oct:days,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m.oct:days+cohort.factor,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m.oct:days+cohort,
			weight~sex*R_POTEMP_25m+sex:days+R_POTEMP_25m.oct:days,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort+stratification,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort+pycnocline_depth,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort+R_POTEMP_25m,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort+R_POTEMP_75m,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort+dynamic_height_0_500m,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort+R_O2_25m,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort+R_O2_75m,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort+R_NO3_25m,
			weight~sex*R_SIGMA_75m+sex:days+R_SIGMA_75m.oct:days+cohort+R_NO3_75m)
} else {
	zcweights.environ.diet=merge(zcweights,data.frame(cohort=1975:lastyear,SST=JunetoSeptAnomalies[4:numyears],SST1=OcttoFebAnomalies[4:numyears],
					SST2=JunetoFebAnomalies[4:numyears],MEI=LaggedMEIJunetoSept[-1],MEI1=LaggedMEIOcttoFeb[-1],
					MEI2=LaggedMEIJunetoFeb[-1],UWI33=UWImeansJunetoSept[1,-(1:6)],UWI36=UWImeansJunetoSept[2,-(1:6)],
					UWI331=UWImeansOcttoFeb[1,-(1:6)],UWI361=UWImeansOcttoFeb[2,-(1:6)],UWI332=UWImeansJunetoFeb[1,-(1:6)],
					UWI362=UWImeansJunetoFeb[2,-(1:6)]))
	fixed.f=list(weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort.factor,
			weight~sex*SST+sex:days+SST1:days+sex:SST1:days+cohort,
			weight~sex*SST+sex:days+SST1:days+sex:SST1:days,
			weight~sex*SST+sex:days+SST1:days+cohort.factor,
			weight~sex*SST+sex:days+SST1:days+cohort,
			weight~sex*SST+sex:days+SST1:days,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days+cohort.factor,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days+cohort,
			weight~sex*UWI36+sex:days+UWI361:days+sex:UWI361:days,
			weight~sex*UWI36+sex:days+UWI361:days+cohort.factor,
			weight~sex*UWI36+sex:days+UWI361:days+cohort,
			weight~sex*UWI36+sex:days+UWI361:days,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days+cohort.factor,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days+cohort,
			weight~sex*MEI+sex:days+MEI1:days+sex:MEI1:days,
			weight~sex*MEI+sex:days+MEI1:days+cohort.factor,
			weight~sex*MEI+sex:days+MEI1:days+cohort,
			weight~sex*MEI+sex:days)
}

zcweights.environ.diet$cohort.factor=factor(ifelse(zcweights.environ.diet$cohort<1990,0,1),labels=c("<=1989",">=1990"))
zcweights.environ.diet$cohort=zcweights.environ.diet$cohort-min(zcweights.environ.diet$cohort)

data(fo)
data(scats)

sardine=get_fo(scats,fo,c("SARSAG"))
anchovy=get_fo(scats,fo,c("ENGMOR"))
sa=get_fo(scats,fo,c("SARSAG","ENGMOR"))
rockfish=get_fo(scats,fo,c("SEBSPP","SEBSP1","SEBSP2"))
hake=get_fo(scats,fo,c("MERPRO"))
squid=get_fo(scats,fo,c("LOLOPA"))
mackerel=get_fo(scats,fo,c("SCOJAP","SCOSPP","SCOMAR","TRAYSYM"))
lo.cal=get_fo(scats,fo,c("SEBSPP","SEBSP1","SEBSP2","MERPRO","LOLOPA"))
hi.cal=get_fo(scats,fo,c("SARSAG","ENGMOR","SCOJAP","SCOSPP","SCOMAR","TRAYSYM"))

fo=cbind(Year=as.numeric(names(sardine)),sardine=sardine,anchovy=anchovy,
		rockfish=rockfish,hake=hake,squid=squid,mackerel=mackerel)
# first and second principal components of 5 primary prey species
fo=cbind(fo,predict(prcomp(fo[,-1]))[,1:2])
fo$lo.cal=lo.cal
fo$hi.cal=hi.cal
fo$ratio=fo$hi.cal/fo$lo.cal
fo$sa=sa

if(use.calcofi)
{
	zcweights.environ.diet$Year=zcweights.environ.diet$cohort+1984
} else{
	zcweights.environ.diet$Year=zcweights.environ.diet$cohort+1975
}
ZCWeight.df1=ZCWeight.df[!is.na(ZCWeight.df$female.environ.mean.fall),]
ZCWeight.df1=ZCWeight.df1[rownames(ZCWeight.df1)%in%rownames(fo),]
zcweights.environ.diet=zcweights.environ.diet[zcweights.environ.diet$Year%in%fo$Year,]
zcweights.environ.diet=merge(zcweights.environ.diet,fo,by="Year")
#
# Next evaluate sequence of fixed-effect models with the chosen random effect model
#
random.f=list(res.environ$best.r)

fixed.f.diet=sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ PC1")))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ squid "))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ rockfish "))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ squid + rockfish "))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ squid + sa "))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ sardine "))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ sa"))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ anchovy"))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ lo.cal"))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ hi.cal"))))
fixed.f.diet=c(fixed.f.diet,sapply(fixed.f,function(x) as.formula(paste("weight~",as.character(x)[3],"+ ratio"))))

res.environ.diet=fitmixed(fixed.f.diet,random.f,data=zcweights.environ.diet) 

#  fit best environment-diet growth model with REML
zc.weight.model.environ.diet=lme(res.environ.diet$best.f,random=res.environ.diet$best.r,data=zcweights.environ.diet,control=lmeControl(opt="optim"),method="REML")
summary(zc.weight.model.environ.diet)

pp=zcweights.environ.diet
pp$days=0
pp1=predict(zc.weight.model.environ.diet,newdata=pp)
pp0=predict(zc.weight.model.environ.diet,newdata=pp,level=0)

pp0=tapply(as.vector(pp0),list(zcweights.environ.diet$cohort,zcweights.environ.diet$sex),mean)
pp1=tapply(as.vector(pp1),list(zcweights.environ.diet$cohort,zcweights.environ.diet$sex),mean)

female.averages=data.frame(fit=pp1[,1])
male.averages=data.frame(fit=pp1[,2])

expected.female.averages=data.frame(fit=pp0[,1])
expected.male.averages=data.frame(fit=pp0[,2])

add=0
if(use.calcofi)add=9
ZCWeight.df1=ZCWeight.df[as.numeric(rownames(ZCWeight.df))%in%(as.numeric(rownames(pp1))+1975+add),]
ZCWeight.df1$female.environ.diet.mean.fall=female.averages$fit
ZCWeight.df1$female.environ.diet.mean.fall.se=female.averages$se
ZCWeight.df1$male.environ.diet.mean.fall=male.averages$fit
ZCWeight.df1$male.environ.diet.mean.fall.se=male.averages$se
ZCWeight.df1$Year=as.numeric(rownames(ZCWeight.df1))


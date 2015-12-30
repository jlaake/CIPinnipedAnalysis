#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
require(CIPinnipedAnalysis)
sdir=system.file(package="CIPinnipedAnalysis")
if(!exists("nboot"))nboot=100
####################################
source(file.path(sdir,"CU_Weight_Adjustment_Model.r"))
source(file.path(sdir,"CU_Weight_Environment_Model.r"))

CUWeight.df=cbind(Year=as.numeric(rownames(CUWeight.df)),CUWeight.df)

store_weights(CUWeight.df,species="CU",fdir=fdir)


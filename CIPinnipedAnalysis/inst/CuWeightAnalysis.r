#if fdir="" it looks for data files in the CalcurData package directory of your R library.
#if fdir=NULL it looks in databases.txt in CalcurData package directory to get the database location
#if fdir is anything else it uses the value of fdir as the directory for database.  
#The scripts check for the value of fdir and if it exists the script will not change the value; otherwise it sets it to NULL
if(!exists("fdir"))fdir=NULL
require(CIPinnipedAnalysis)
sdir=system.file(package="CIPinnipedAnalysis")
if(!exists("nboot"))nboot=100
####################################
# Set this value; be aware that all of the environmental data has to be entered through Feb of lastyear+1 
# for the growth script to work properly; script CreateAnaomalies.r uses lastyear;
# change value of lastyear prior to running scripts and it will override the
# setting below.
if(!exists("lastyear"))lastyear=2014
####################################
source(file.path(sdir,"CU_Weight_Adjustment_Model.r"))
source(file.path(sdir,"CU_Weight_Environment_Model.r"))

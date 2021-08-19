# The following script will simplify the simulations process
# moving forward, by automatically generating all necessary
# scripts to run simulations and dynamically adjusting 
# variables in them.

# Please enter the location of the folder from which you are
# running this script (which should also contain the xDatasets
# and xInputFiles subfolders)
ScriptFolder <- '/Volumes/T3/mplusanalyzer'

# Specify the TARGET folder where you would like to generate
# the analysis setup. This does not need to be where you
# eventually run it!
TargetFolderRoot <- '/Volumes/T3/local'

# Specify the path to the location where you will be
# RUNNING the scripts. If doing it on a cluster, then
# this would be the path to where the folder will be in
# your cluster.
#RunningFolderRoot <- '/scratch/calvacfa'
RunningFolderRoot <- TargetFolderRoot

# Specify sample ID
# This is VERY important since it will define which functions will
# be used. It ~MUST~ be one of the following:
# - 'NESARCw1' (NESARC wave 1)
# - 'NESARCw2' (NESARC wave 2)
# - 'NCS' (National Comorbidity Survey)
# - 'NCSR' (National Comorbidity Survey - Replication)
# This will also serve as a prefix  for results CSVs and graphs.
SampleID <- 'NESARCw1'

# Specify whether you'll do a CFA only, CFA with EFA, or EFA only
# Valid options: "CFA", "EFA", "CFA+EFA"
AnalysisType <- "CFA"

# Specify model to examine
ModelID <- '2Fvf'

# Save difftest .dat file?
SaveDiffTest <- TRUE

# Do difftest?
DoDiffTest <- FALSE

# Specify which ways you would like to drop zero cases
# (i.e., what guiding principle for dropping them)
DROP_OVERALL_ZERO <- TRUE # In NESARC, this is countdx
DROP_FACTORS <- FALSE # By the sum of all diagnoses for a given factor
DROP_BY_DX <- FALSE # By each specific diagnostic variable

# Specify number of simulations to run per condition
# I like to use 1 for testing, and 500 for production.
NumSimulationsPerCondition <- 500

# Do we want to run the external covariates step?
RUN_EXT_COVARS <- FALSE # Valid inputs: TRUE or FALSE # IGNORE THIS
# Specify demographics columns to use for ext analyses
Demos <- c("Demo1") # leave this as is if not running that step

############################################################
#### DON'T CHANGE ANYTHING BELOW HERE!! ####################
############################################################

if(SaveDiffTest & DoDiffTest) {
  print("YOU CAN'T HAVE SAVE DIFFTEST AND DO DIFFTEST BOTH AT SAME TIME!")
  stop()
}


# Out folder name is automatically computed now
outFolderName <- sprintf("%s_%s_%d"
                         ,SampleID,ModelID
                         ,NumSimulationsPerCondition)

# And thus, so are RunningFolder and TargetFolder
RunningFolder <- file.path(RunningFolderRoot,outFolderName)
TargetFolder <- file.path(TargetFolderRoot,outFolderName)

# Driver Script Template filename
Filename_DriverScriptTmp <- 'DriverScript_Template_MM.R'
# Scraping functions filename
Filename_ScrapingFunctions <- 'ScrapingFunctions.R'
# Implemented list filename
Filename_Globals <- 'Globals.R'
# ExtModel Class File --- IGNORE ALL THIS
Filename_EMclass <- 'ExtModelClass.R'
Filename_EMfuns <- 'ExtModelFunctions.R'
Filename_InpFuns <- "InputFileFunctions.R"
Filename_EMtemplate <- 'RunEV_Template.R'

library(tidyverse)
library(readxl)
library(arrangements)
source(file.path(ScriptFolder,Filename_Globals))
source(file.path(ScriptFolder,Filename_ScrapingFunctions))

# Validate SampleID!
validSamples <- VALID_SAMPLES
if (!(SampleID %in% validSamples)) {
  print(sprintf('ERROR: Invalid sample ID (you stated %s)'
                ,SampleID))
  stop()
}

# Compute VALID_INPUTS global
ModSpecFilename <- sprintf("%s.xlsx",SampleID)
PathXLSX <- file.path(ScriptFolder,'xModelSpecs',ModSpecFilename)
globals <- ComputeGlobals(SampleID,PathXLSX)
VALID_INPUTS <- globals[['VALID_INPUTS']]
MM_INPFILES <- globals[['MM_INPFILES']]
#SM_INPFILES <- globals[['SM_INPFILES']]

# Validate ModelID!
validModels <- VALID_INPUTS[[SampleID]]
if (!(ModelID %in% validModels)) {
  print(sprintf('ERROR: Invalid model ID (you stated %s)'
                ,ModelID))
  stop()
}

# Specify the name of the dataset file to be used:
# (this file should be located in ScriptFolder/xDatasets)
DatasetFile <- DATASET_NAMES[[SampleID]]


# Note that we can't really od the "drop factors" thing if we only have
# one factor!
num_factors <- length(names(MODEL_SPECS[[SampleID]][[ModelID]][['facs']]))
if (num_factors == 1) {
  DROP_FACTORS <- FALSE
}



# Ok, now validate if there's a valid definition specified.
if (!(DROP_OVERALL_ZERO | DROP_FACTORS | DROP_BY_DX)) {
  print('ERROR: You must specify at least one valid definiton for dropping 0 cases.')
  if (num_factors == 1) {
    print ('(Note: you cannot use the DROP_FACTORS option if only one factor was specified)')
  }
  stop()
}

# Read model structure
ModelTable <- read_excel(file.path(ScriptFolder,'xModelSpecs',ModSpecFilename)
                         , sheet="CSV_Structure")
ModelStructure <- ComputeModelStructCW(SampleID,ModelTable,KEY_VARS)
MODEL_SPECS[[SampleID]] <- ModelStructure[['ModelStruct']]
CROSSWALKS[[SampleID]] <- ModelStructure[['Crosswalk']]

# Specify the name of the Mplus suffix file to be used:
# (this file should be located in ScriptFolder/xInputFiles)
if (DoDiffTest) {
  SuffixFileMeasurementModel <- gsub("\\.inp","_Diff.inp",MM_INPFILES[[SampleID]][[ModelID]])
  SuffixFileMeasurementModel_Bak <- MM_INPFILES[[SampleID]][[ModelID]]
} else {
  SuffixFileMeasurementModel <- MM_INPFILES[[SampleID]][[ModelID]]
}

# This is the file to be used for step 4, if performed
#SuffixFileStructuralModel <- SM_INPFILES[[SampleID]][[ModelID]]

# TODO: fix this shit!!!
# Intelligently determine the splitting variables.
SplitVariables <- c()
if (DROP_OVERALL_ZERO) {
  SplitVariables <- c(SplitVariables,
                      FacSplitNames(SampleID,ModelID,MODEL_SPECS
                                    ,KEY_VARS,'overall'))
}
if (DROP_FACTORS & !(str_detect(ModelID,"1F"))) {
  SplitVariables <- c(SplitVariables,
                      FacSplitNames(SampleID,ModelID,MODEL_SPECS
                                    ,KEY_VARS,'facs'))
}
if (DROP_BY_DX) {
  SplitVariables <- c(SplitVariables,
                      FacSplitNames(SampleID,ModelID,MODEL_SPECS
                                    ,KEY_VARS,'dx'))
}

# Ok, so now we will prepare this for you.

# Navigate to directory
setwd(ScriptFolder)

# First, we will read the template script file.
con <- file(Filename_DriverScriptTmp,open='r')
ScriptTemplate <- readLines(con)
close(con)

# Replace the general variables that will not change
ScriptTmpGen <- gsub('GSUB_SAMPLEID'
                     ,SampleID
                     ,ScriptTemplate
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_MODELID'
                     ,ModelID
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_WDPATH'
                     ,RunningFolder
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_MM_SUFFIXFILE'
                     ,SuffixFileMeasurementModel
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
if (DoDiffTest) {
  ScriptTmpGen <- gsub('GSUB_MMBAK_SUFFIXFILE'
                       ,SuffixFileMeasurementModel_Bak
                       ,ScriptTmpGen
                       ,ignore.case = FALSE)
}
ScriptTmpGen <- gsub('GSUB_RUNTYPE'
                     ,AnalysisType
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_DATAFILE'
                     ,DatasetFile
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_NUMSIMS'
                     ,NumSimulationsPerCondition
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_DOEXTVARS'
                     ,ifelse(RUN_EXT_COVARS,"TRUE","FALSE")
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_DROPFACTORS'
                     ,ifelse(DROP_FACTORS,
                             ifelse(str_detect(ModelID,'1F')
                                    ,'FALSE','TRUE')
                             ,'FALSE')
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_MODSPECSFILE'
                     ,ModSpecFilename
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_DEMOS'
                     ,cv_to_str(Demos)
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_EVTEMPLATEFILENAME'
                     ,Filename_EMtemplate
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_SAVEDIFFTEST'
                     ,ifelse(SaveDiffTest,"TRUE","FALSE")
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)
ScriptTmpGen <- gsub('GSUB_PERFORMDIFFTEST'
                     ,ifelse(DoDiffTest,"TRUE","FALSE")
                     ,ScriptTmpGen
                     ,ignore.case = FALSE)


# Generate a list of scripts that will need to be generated
if (AnalysisType == "EFA") {
  ScriptList <- sprintf('%s-%s-EFA.R',SampleID,SplitVariables)
} else {
  ScriptList <- sprintf('%s-%s.R',SampleID,SplitVariables)
}

# Generate target folder
dir.create(TargetFolder)

# Copy dataset and Mplus input suffix file to target folder
tmp.filestocopy <- c(file.path(ScriptFolder,'xInputFiles',SuffixFileMeasurementModel)
                     ,file.path(ScriptFolder,'xDatasets',DatasetFile)
                     ,file.path(ScriptFolder,Filename_InpFuns))
if (AnalysisType %in% c("EFA","CFA+EFA")) {
  if (DoDiffTest) {
    efa.filename <- gsub("_MM_Diff","_EFA",SuffixFileMeasurementModel)
  } else {
    efa.filename <- gsub("_MM","_EFA",SuffixFileMeasurementModel)
  }
  tmp.filestocopy <- c(tmp.filestocopy
                       ,file.path(ScriptFolder,'xInputFiles',efa.filename))
}
if (RUN_EXT_COVARS) {
  tmp.filestocopy <- c(tmp.filestocopy
                       ,file.path(ScriptFolder,Filename_EMclass)
                       ,file.path(ScriptFolder,Filename_EMfuns)
                       ,file.path(ScriptFolder,Filename_EMtemplate)
  )
}
if (DoDiffTest) {
  # copy the file without difftest as backup
  tmp.filestocopy <- c(tmp.filestocopy,file.path(ScriptFolder,'xInputFiles',SuffixFileMeasurementModel_Bak))
}

file.copy(tmp.filestocopy,TargetFolder)

# copy model specs file
  file.copy(file.path(ScriptFolder,'xModelSpecs',ModSpecFilename)
            ,TargetFolder)

# Copy scraping functions and globals auxiliary file to be used
file.copy(c(file.path(ScriptFolder,Filename_ScrapingFunctions)
            ,file.path(ScriptFolder,Filename_Globals))
          ,TargetFolder)

# For each of our 'zero variables' we want to run simulations
# for, generate a separate script
for (i in 1:length(SplitVariables)) {
  # Fill it out
  myScript <- gsub('GSUB_ZEROVARIABLE'
                   ,SplitVariables[i]
                   ,ScriptTmpGen
                   ,ignore.case = FALSE)
  # Save the file
  print(sprintf('Saving script file %s...',ScriptList[i]))
  myCon <- file(file.path(TargetFolder,ScriptList[i]))
  writeLines(myScript,con = myCon)
  close(myCon)
  # prep
  dir.create(file.path(TargetFolder,SplitVariables[i]))
}

# Fantastic. Now, let's finish this up by making a script
# to automate the process in case we want to batch run it :)
LogList <- gsub('.R','.log',ScriptList)
ScriptLines <- sprintf('Rscript %s |& tee %s'
                       ,ScriptList,LogList)
ScriptLines <- c("#!/bin/tcsh","",ScriptLines,"")
sCon <- file(file.path(TargetFolder,'RunAll.sh'))
writeLines(ScriptLines,con = sCon)
close(sCon)

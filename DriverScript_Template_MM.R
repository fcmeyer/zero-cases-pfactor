###########################################
# Settings                                #
###########################################\

# Set seed for replication purposes
SEED <- 031289

# Specify sample ID
# This is important since it will define which functions will
# be used.
SAMPLE_ID <- 'GSUB_SAMPLEID'

# Specify model ID
# This is important since it will define which functions will
# be used.
MODEL_ID <- 'GSUB_MODELID'

# Specify name of measurement model suffix file 
# (present in the directory specified above)
MM_SUFFIX_FILE <- 'GSUB_MM_SUFFIXFILE'

# Indicate number of simulations to run for each condition
N_SIMS_PER_COND <- GSUB_NUMSIMS
# Indicate the "step" by which percentages of dropped
# cases will grow.
STEP <- 0.10
# Indicate the minimum and maximum values for percentages of
# zero cases retained
START_VALUE <- 0
END_VALUE <- 1.00

# Specify root working directory here:
WD_PATH <- 'GSUB_WDPATH'

# Specify zero control variable
ZERO_VAR <- 'GSUB_ZEROVARIABLE'

# Specify path to dataset file
DATA_FILE <- 'GSUB_DATAFILE'

# Specify path to model specs file (if applicable)
MODSPEC_FILE <- 'GSUB_MODSPECSFILE'

# Compute drop factors variables?
DROP_FACTORS <- GSUB_DROPFACTORS

# Specify if we're running EFA or not
RUN_TYPE <- "GSUB_RUNTYPE"

# Specify if Step 4 is completed
RUN_EXT_COVARS <- GSUB_DOEXTVARS

# Specify if we want to save difftest
SAVE_DIFFTEST <- GSUB_SAVEDIFFTEST

# Specify if we want to perform difftest
PERFORM_DIFFTEST <- GSUB_PERFORMDIFFTEST

# ignore this, it's for another project
if (RUN_EXT_COVARS) {
  DEMOS <- GSUB_DEMOS
  Filename_EVTmp <- "GSUB_EVTEMPLATEFILENAME"
}

if (PERFORM_DIFFTEST) {
  BAK_MM_SUFFIX_FILE <- "GSUB_MMBAK_SUFFIXFILE"
}

###########################################
# Prepare to run scripts! (ALWAYS RUN!)   #
###########################################

# Load required libraries
library(tidyverse)
library(dplyr)
library(MplusAutomation)
library(texreg)
library(parallel)
library(doParallel)
library(foreach)
library(arrangements)
library(readxl)
library(R6)
library(hablar)

# Change working directory to our folder.
setwd(WD_PATH)
# Set our replication seed.
set.seed(SEED)
# Load the Scraping Functions!
source(file.path(WD_PATH,'ScrapingFunctions.R'))
# Load the Globals!
source(file.path(WD_PATH,'Globals.R'))
# load input file funs
source(file.path(WD_PATH,'InputFileFunctions.R'))
# Load the Ext Model Stuff (if needed)
# IGNORE THIS -- for another project
if(RUN_EXT_COVARS) {
  source(file.path(WD_PATH,'ExtModelClass.R'))
  source(file.path(WD_PATH,'ExtModelFunctions.R'))
}


# This is assumed
if (RUN_TYPE %in% c("CFA+EFA","EFA")) {
  if (PERFORM_DIFFTEST) {
    EFA_SUFFIX_FILE <- gsub("_MM_Diff","_EFA",MM_SUFFIX_FILE)
  } else {
    EFA_SUFFIX_FILE <- gsub("_MM","_EFA",MM_SUFFIX_FILE)
  }
}

# adding this for the free variance files
SCRAPE_FAC_VARIANCE <- FALSE
if (str_detect(MODEL_ID,"vf")) {
  SCRAPE_FAC_VARIANCE <- TRUE
}
  

###########################################
# Stage 1: Data Prep                      #
###########################################
# This stage will prepare the input files #
# that will be used by Mplus. Essentially #
# doing everything prior to running the   #
# models themselves.                      #
###########################################

# Read appropriate data file
data <- read_csv(file.path(WD_PATH,DATA_FILE), col_names=TRUE ,na='.')

# We'll fix SIM_PREFIX to the name of ZERO_VAR for now.
SIM_PREFIX <- ZERO_VAR

# Compute VALID_INPUTS global
globals <- ComputeGlobals(SAMPLE_ID,file.path(WD_PATH,MODSPEC_FILE))
VALID_INPUTS <- globals[['VALID_INPUTS']]
MM_INPFILES <- globals[['MM_INPFILES']]

# Enhance or calculate globals as needed.
ModelTable <- read_excel(file.path(WD_PATH,MODSPEC_FILE)
                         ,sheet="CSV_Structure")
ncs_struct <- ComputeModelStructCW(SAMPLE_ID,ModelTable,KEY_VARS)
MODEL_SPECS[[SAMPLE_ID]] <- ncs_struct[['ModelStruct']]
CROSSWALKS[[SAMPLE_ID]] <- ncs_struct[['Crosswalk']]

# If we're doing the external covariate analyses, read that too!
if (RUN_EXT_COVARS) {
  ExtVarTable <- read_excel(file.path(WD_PATH,MODSPEC_FILE)
                           ,sheet="CSV_ExternalVars")
}

# Adjust start value if it's a single diagnosis
if(ZERO_VAR %in% CROSSWALKS[[SAMPLE_ID]]$oglabs) {
  START_VALUE <- STEP
}

# In the following block, we will set up a dataframe to track 
# and label each simulation that will be run.
# -----------------------------------------------------------
# cond: one row for each simulation, specify condition as 
#       fraction (eg, .1 .1 .1 .2 .2 .2 ... .9 .9 .9) if
#       N_SIMS_PER_COND = 3

# Compute a list of conditions (based on specified settings)
conditions <- seq(START_VALUE,END_VALUE,STEP)

cond <- rep(conditions
            ,times=(rep(N_SIMS_PER_COND,length(conditions))))

# cond_pct: same as above, but in percentage form
cond_pct <- cond*100
# cond_ID: this is going to be the ID within each condition
cond_ID <- rep(1:N_SIMS_PER_COND
               ,times=length(conditions))
# this is the overall id, which concatenates cond with cond_ID
sim_ID <- sprintf("%03.f_%03.f"
                  ,cond_pct
                  ,cond_ID)
# the names of each input file
inp_file <- sprintf("Sim_%s-%s.inp",SIM_PREFIX,sim_ID)
# the names of each corresponding dat file
dat_file <- sprintf("Sim_%s-%s.dat",SIM_PREFIX,sim_ID)
# the names of each corresponding out file
out_file <- sprintf("Sim_%s-%s.out",SIM_PREFIX,sim_ID)

# compile all the above into a new dataframe!
Simulations <- data.frame(cond_pct,cond_ID,sim_ID
                          ,inp_file,dat_file,out_file
                          ,stringsAsFactors = FALSE)

# For EFA only:
if (RUN_TYPE %in% c("EFA","CFA+EFA")) {
  Simulations$efa.inp_file <- sprintf("Sim_%s-%s_EFA.inp",SIM_PREFIX,sim_ID)
  Simulations$efa.out_file <- sprintf("Sim_%s-%s_EFA.out",SIM_PREFIX,sim_ID)
  Simulations$efa.gh5_file <- sprintf("Sim_%s-%s_EFA.gh5",SIM_PREFIX,sim_ID)
}

# We only need ONE simulation with 0% of asymptomatic cases
# since it would always turn out the same!
if (0 %in% conditions) {
  Simulations <- Simulations %>% group_by(cond_pct) %>% filter(cond_pct > 0 | cond_ID == 1)
  Simulations <- as.data.frame(Simulations)
}
# We only need ONE simulation with 100% of asymptomatic cases
# since it would always turn out the same! (i.e., full sample)
if (1 %in% conditions) {
  Simulations <- Simulations %>% group_by(cond_pct) %>% filter(cond_pct < 100 | cond_ID == 1)
  Simulations <- as.data.frame(Simulations)
}


# Perform data cleanup operations (remove useless columns)
keepVars <- MODEL_SPECS[[SAMPLE_ID]][[MODEL_ID]][['mmVars']]
# TODO: ensure I adjust this correctly for when not just doing CFA.
if (RUN_EXT_COVARS) {
  # We also want to keep the relevant covariates in our dataset.
  keepVars <- c(keepVars,ExtVarTable$Label)
}
data <- data %>% select(all_of(keepVars))

# Create var that's sum of all diagnoses
data <- computeSumDx(data,SAMPLE_ID,MODEL_ID
                     ,MODEL_SPECS,KEY_VARS
                     ,VALID_INPUTS)

# Create new variable that computes sum of internalizing
# and externalizing diagnoses
if (DROP_FACTORS) {
  data <- computeFacSplit(data,SAMPLE_ID,MODEL_ID
                          ,MODEL_SPECS,CROSSWALKS
                          ,VALID_INPUTS)
}

# Grab our filter variable!
FilterVar <- data[,ZERO_VAR]

# Create two data frames:
# 0. All zero cases
mat_0 <- as.data.frame(data[which(FilterVar == 0),])
# 1. All cases with at least one dx 
mat_1 <- as.data.frame(data[which(FilterVar > 0),])

# Compute the total number of zero cases in sample
num_0dx <- dim(mat_0)[1]

# Great! Now, let's generate ALL these files to prepare
# in advance for running our simulations. Yee haw!

# First, read in the suffix file, which is same for all
# cases.
con <- file(MM_SUFFIX_FILE,open='r')
MplusSpecs <- readLines(con)
close(con)

# Read in backup input file in case we don't have a difftest file
if (PERFORM_DIFFTEST) {
  conBk <- file(BAK_MM_SUFFIX_FILE,open='r')
  MplusSpecsBK <- readLines(conBk)
  close(conBk)
}

# And read the EFA suffix file, if applicable
if (RUN_TYPE %in% c("EFA","CFA+EFA")) {
  conx <- file(EFA_SUFFIX_FILE,open='r')
  MplusSpecs.EFA <- readLines(conx)
  close(conx)
}

# Ok, now we want to create a new sub-directory
# for this specific set of simulations, based on SIM_PREFIX.
SimDirPath <- file.path(getwd(),SIM_PREFIX)
dir.create(SimDirPath)
setwd(SimDirPath)

# For saving Mplus file headers (only used in Ext analysis)
MplusHeaderList <- list()

# Now, let's generate input and data files for each of the
# simulations we will be running!
for (i in 1:length(Simulations$sim_ID)) {
  # Get me a condition fraction for multiplying!
  i_cond <- Simulations$cond_pct[i] / 100
  # Number of observations to keep
  samp_int <- as.integer(round(i_cond*num_0dx))
  # Don't do any binding if we don't keep any zero cases!
  if (i_cond == 0) {
    full_sample <- mat_1
  } else if (i_cond == 1) {
    full_sample <- as.data.frame(data)
  } else {
    # Sample my 0 diagnosis cases (withOUT replacement!)
    zeroDx_sample <- sample_n(mat_0, samp_int, replace = FALSE)
    # Merge them with my 1+ diagnosis cases!
    full_sample <- rbind(zeroDx_sample,mat_1)
  }
  if (RUN_TYPE == "EFA") {
    MplusHeader <- prepareMplusData(full_sample
                                    ,Simulations$dat_file[i]
                                    ,inpfile=FALSE)
  } else {
    # Generate an Mplus input and simulation-specific data file 
    MplusHeader <- prepareMplusData(full_sample
                                    ,Simulations$dat_file[i]
                                    ,inpfile=FALSE)
  }
  
  if (RUN_TYPE %in% c("CFA+EFA","CFA")) {
    
    # Write out an input file that concatenates the prefix
    # from above and our suffix file!
    con2 <- file(Simulations$inp_file[i],open='w')
    write.tmp <- unlist(c(MplusHeader,MplusSpecs))
    
    # if we're performing a difftest, ensure we fill that in with our
    # files!
    if (PERFORM_DIFFTEST) {
      diff_filename <- gsub("\\.inp","-diff.dat",Simulations$inp_file[i])
      # ensure file exists, if not stop.
      if(file.exists(file.path(WD_PATH,ZERO_VAR,diff_filename))) {
        write.tmp <- gsub("GSUB_DIFFFILENAME",diff_filename,write.tmp)
      } else {
        # run the model without the difftest
        write.tmp <- unlist(c(MplusHeader,MplusSpecsBK))
      }
    }
    # If we're running a second step and looking to use factor
    # scores, make sure we tack that onto the model!
    #if (RUN_EXT_COVARS & (E.METHOD == "BOTH" | E.METHOD == "FSE")) {
    # By default, now save FSEs and difftest stuff. More outputs, but whatever.
    SaveSuffix <- L.SaveFS(SIM_PREFIX,Simulations$sim_ID[i])
    write.tmp <- unlist(c(write.tmp,SaveSuffix))
    
    if (SAVE_DIFFTEST) {
      diff_filename <- gsub("\\.inp","-diff.dat",Simulations$inp_file[i])
      diff.line <- sprintf("DIFFTEST IS %s;",diff_filename)
      write.tmp <- c(write.tmp,"",diff.line,"")
    }
    
    #}
    writeLines(write.tmp,con = con2)
    close(con2)
  }
  
  if (RUN_TYPE %in% c("CFA+EFA","EFA")) {
    # Write out an input file for the EFA portion, concatenating
    # the prefix  from above and our suffix file!
    conb <- file(Simulations$efa.inp_file[i],open='w')
    write.tmp <- unlist(c(MplusHeader,MplusSpecs.EFA))
    writeLines(write.tmp,con = conb)
    close(conb)
  }
  
  # If we are doing external analyses, we should also save the headers. They
  # will come in super handy later!
  # IGNORE - FOR ANOTHER PROJECT
  if (RUN_EXT_COVARS) {
    MplusHeaderList[[Simulations$sim_ID[i]]] <- MplusHeader
  }
  
}

# Let's save our simulations table too, just for future-
# proofing purposes.
write.csv(Simulations
          ,sprintf('%s-%s_SimulationsList.csv',SAMPLE_ID,SIM_PREFIX)
          ,row.names = FALSE)

# We now have all of our input and data files! It is time to
# start running stuff. 

###########################################
# Stage 2: Running simulations            #
###########################################

# Run all the Mplus models at once! Yee haw!
# WARNING: this will take a while!!!!
# Now in parallel! YAYYYYY

if (RUN_TYPE == "EFA") {
  file_list <- list.files(path=SimDirPath
                          ,pattern = '.*EFA.inp')
} else {
  file_list <- list.files(path=SimDirPath
                          ,pattern = '.*inp')
}

numCores <- detectCores()
registerDoParallel(numCores)
rv <- foreach(i=1:length(file_list)) %dopar% {
  runModels(target = file_list[i]
            ,showOutput=FALSE
            ,local_tmpdir = TRUE
            ,logFile = NULL)
}

###########################################
# Stage 3: Scraping simulation data       #
###########################################

# Number of cases?
n_sims <- length(Simulations$sim_ID)

# Following, we will create a results dataframe to contain
# all the data we want to scrape out of our models.
SimResults <- createContainerDF(Simulations,SAMPLE_ID,MODEL_ID
                                ,VALID_INPUTS,MODEL_SPECS,RUN_TYPE,PERFORM_DIFFTEST)

# Get me the structure of the model
facs <- MODEL_SPECS[[SAMPLE_ID]][[MODEL_ID]][['facs']]
phis <- MODEL_SPECS[[SAMPLE_ID]][[MODEL_ID]]$phis

# Now, we will scrape the results of the CFA. Note that here, we are
# parallelizing the task by splitting the dataframe into
# rows, scraping each row (simulation) separately, and then
# merging them all together at the end.
if (RUN_TYPE %in% c("CFA","CFA+EFA")) {
  SimResults <- foreach(i = 1:nrow(SimResults)
                        , .combine=rbind) %dopar% {
                          tmp <- scrapeSim( SimResults[i,]
                                            ,SIM_PREFIX
                                            ,SAMPLE_ID
                                            ,MODEL_ID
                                            ,VALID_INPUTS
                                            ,CROSSWALKS
                                            ,MODEL_SPECS
                                            ,PERFORM_DIFFTEST
                          )
                          tmp
                        }
  
  # Oops! forgot to include the condition and ids (to make it easier
  # for plots). Let's get that back!
  tmp <- Simulations[,c('cond_pct','cond_ID','sim_ID')]
  names(tmp)[3] <- 'sids'
  Sims <- inner_join(tmp
                     ,SimResults
                     ,by = 'sids')
  
  # One last thing: we're going to also want to compute
  # H reliability index for each factor we estimated.
  Sims <- compute_Hindex(Sims,facs)
  
  # Compute omega
  Sims$omega <- compute_omega(Sims,facs)
  Sims <- compute_omegaH_relomega(Sims,facs)
  
  # If bifactor, compute the other ones.
  if (str_detect(MODEL_ID,"Bi")) {
    Sims$ecv <- compute_ecv(Sims,facs)
    Sims <- compute_omegaHS(Sims,facs)
    puc <- compute_puc(facs)
    Sims$puc <- puc
  }
  
  # Ok, now save that again.
  write.csv(Sims
            ,sprintf('%s-%s_SimulationsResults.csv'
                     ,SAMPLE_ID,SIM_PREFIX)
            ,row.names = FALSE)
  
} else if (RUN_TYPE == "EFA") {
  Sims <- read_csv(sprintf('%s-%s_SimulationsResults.csv'
                           ,SAMPLE_ID,SIM_PREFIX),col_names=TRUE)
}

###########################################
# Stage 4: EFA                            #
########################################### 

if (RUN_TYPE %in% c("CFA+EFA","EFA")) {
  
  SimResults.EFA <- createContainerDF(Simulations,SAMPLE_ID,MODEL_ID
                                      ,VALID_INPUTS,MODEL_SPECS,"EFA")
  
  SimResults.EFA <- foreach(i = 1:nrow(SimResults.EFA)
                            , .combine=rbind) %dopar% {
                              tmp <- scrapeSim.EFA( SimResults.EFA[i,]
                                                    ,SIM_PREFIX
                              )
                              tmp
                            }
  
  # Oops! forgot to include the condition and ids (to make it easier
  # for plots). Let's get that back!
  tmp <- Simulations[,c('cond_pct','cond_ID','sim_ID')]
  names(tmp)[3] <- 'sids'
  Sims.EFA <- inner_join(tmp
                         ,SimResults.EFA
                         ,by = 'sids')
  # Move the error column to last, it's annoying!
  Sims.EFA <- Sims.EFA %>% select(-errs.EFA,errs.EFA)
  
  # Ok, now save the EFA only results.
  write.csv(Sims.EFA
            ,sprintf('%s-%s_EFA_SimulationsResults.csv'
                     ,SAMPLE_ID,SIM_PREFIX)
            ,row.names = FALSE)
  
  # And now, let's also combine the EFA stuff with the CFA stuff into 
  # one table, and save that!
  Sims.Combined <- inner_join(Sims,Sims.EFA,by=c("cond_pct","cond_ID","sids"))
  write.csv(Sims.Combined
            ,sprintf('%s-%s_CFAEFA_SimulationsResults.csv'
                     ,SAMPLE_ID,SIM_PREFIX)
            ,row.names = FALSE)
}


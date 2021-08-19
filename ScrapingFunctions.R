# This file contains useful functions for scraping
# Mplus outputs.


# Compute globals from model specs ----------------------------------------

# Model specifications for each sample are given in the ModelSpecs file
# for a given sample.
ComputeGlobals <- function(SampleID,PathXLSX) {
  
  InputFileTab <- read_excel(PathXLSX,sheet="CSV_InputFiles")
  
  # Valid inputs for this sample
  ValidInputs <- list()
  ValidInputs[[SampleID]] <- InputFileTab$Model_ID
  
  # Models and corresponding imput files
  # Measurement models and EFA (MM)
  MM.Sample <- InputFileTab %>% 
                 deframe() %>% 
                 as.list()
  MM <- list()
  MM[[SampleID]] <- MM.Sample
  
  # Return all that to the user
  rv <- list(VALID_INPUTS = ValidInputs
             ,MM_INPFILES = MM
             #,SM_INPFILES = SM
             )
  return(rv)
}

# This function will compute the module structure for a complex sample,
# in a way that results in equivalent format to the one you can see was
# manually specified in the Globals file.
# Inptus required are as follows:
#   - SampleID: how you'd like to call the sample
#   - ModelTable: table MUST include Label, FM_Label, and model columns.
# Outputs:
#   - a list containing the Crosswalk and the ModelSpecs for adding to
#     globals!
# IMPORTANT NOTE: library(arrangement) is required for this function!!
ComputeModelStructCW <- function(SampleID,ModelTable,KeyVars) {
  # Remove labels in KEY_VARS to avoid problems!
  for(sample in names(KeyVars)) {
    names(KeyVars[[sample]]) <- NULL
  }
  
  # Compute the crosswalk table for this sample
  Crosswalk <- ModelTable %>% 
                select(Label,FM_Label) %>%
                rename(oglabs = Label
                       ,mylabs = FM_Label)
  # Now, compute the model structure
  # Get the list of models
  ModelList <- ModelTable %>%
                select((!contains("Label")),-"VarType") %>%
                names()
  ModelStruct <- list()
  for (model in ModelList) {
    # Determine if model is bifactor
    isBiFactor <- str_detect(model,"Bi.")
    # Slice our model table for the model in question
    tmpMT <- ModelTable %>% 
              select(Label,FM_Label,all_of(model)) %>% 
              rename(fctr = all_of(model)
                     ,mylab = FM_Label
                     ,oglab = Label) %>%
              drop_na() %>% 
              mutate_at(vars(fctr)
                        ,str_sub,start=1,end=1)
    # Identify factors present
    fac_names <- unique(tmpMT$fctr)
    # Construct facs list
    facs <- list()
    for (fac_name in fac_names) {
      tmpMTfac <- tmpMT %>% filter(fctr == fac_name)
      facs[[fac_name]] <- tmpMTfac$mylab
    }
    # If it's a bifactor model:
    # (1) Need to add GF!
    #     (This is assumed to be ALL variables in the model)
    # (2) No correlations (orthogonal!)
    if (isBiFactor) {
      facs[["G"]] <- tmpMT$mylab
      phis <- NA
    } else { # Not orthogonal; compute correlations!
      if (length(fac_names) == 1) {
        # Only one factor; no correlations necessary
        phis <- NA
      } else {
        # We need to add "phis" statements for all possible
        # correlations.
        cmb <- combinations(fac_names, k=2, replace=FALSE)
        phis <- apply(cmb,1,paste,sep='',collapse='')
      }
    }
    # The last thing we'll do, is add a list of the required
    # variables for modeling the measurement model. 
    # This will help trim the dataset size
    # for each simulation, saving hard drive space.
    mmVars <- c(KeyVars[[SampleID]],tmpMT$oglab)
    
    # Now, get me the "usevars" for this measurement model
    useVars <- tmpMT$oglab
    
    # Create a table with variable types
    varTypes <- ModelTable %>% 
                  select(Label,FM_Label,VarType) %>%
                  rename( mylab = FM_Label
                         ,oglab = Label
                         ,type  = VarType)
    # This only includes the usevars we have though.
    varTypes <- varTypes[varTypes$oglab %in% mmVars,]
    
    ModelStruct[[model]] <- list("facs" = facs
                                 ,"phis" = phis
                                 ,"mmVars" = mmVars
                                 ,"useVars" = useVars
                                 ,"varTypes" = varTypes)
  }
  return(list("Crosswalk" = Crosswalk
              ,"ModelStruct" = ModelStruct))
}

# Function to compute structural model properties
# This will generate a structure that drives external variable tests
ComputeSMstruct <- function(SampleID,PathXLSX) {
  
  SM_VarTab <- read_excel(PathXLSX,sheet="CSV_ExternalVars")

  # Which columns are demographics?
  demo.names <- names(SM_VarTab)[str_detect(names(SM_VarTab),"Demo")]
  # Validate that at least one demographic column is present
  if (length(demo.names) == 0) {
    print("No demographics column (Demo1) present in model. Pls fix and rerun!")
    stop()
  }
  
  # Compile list of lists of demographic variables
  Demos <- list()
  for (demo.name in demo.names) {
    Demos[[demo.name]] <- SM_VarTab$Label[pull(SM_VarTab
                                               ,demo.name)]
  }
  
  # Indicate all external variables specified
  ExtVars <- SM_VarTab$Label[pull(SM_VarTab,"ExtVar")]
  
  # Compile list of all vars that are categorical
  CatVars <- SM_VarTab$Label[pull(SM_VarTab,"isCategorical")]
  
  return(list("Demos" = Demos
              ,"ExtVars" = ExtVars
              ,"Categorical" = CatVars))
}


# Compute overall sum of diagnoses ----------------------------------------

computeSumDx <- function(d # dataset
                            ,SampleID # sample ID
                            ,ModelID # model being run
                            ,ModelSpecs # from global
                            ,KeyVars # from global
                            ,valid.Inputs # from global
) {
  
  # Remove labels in KEY_VARS to avoid problems!
  for(sample in names(KeyVars)) {
    names(KeyVars[[sample]]) <- NULL
  }
  
  rv <- NULL
  valid.Samples <- names(valid.Inputs)
  if (!(SampleID %in% valid.Samples)) {
    print(sprintf('INVALID SAMPLE SPECIFIED: %s! QUITTING.',SampleID))
  } else if (!(ModelID %in% valid.Inputs[[SampleID]])) {
    print(sprintf('INVALID MODEL SPECIFIED FOR SAMPLE %s: %s',SampleID,ModelID))
    print(sprintf('VALID MODELS FOR THIS SAMPLE ARE %s'
                  ,paste(valid.Inputs[[SampleID]],collapse=' ')))
    print('QUITTING.')
  } else {
    # Variable names are formulaic.
    if ((SampleID == 'NESARCw1') | (SampleID == 'NCS')) {
      varNameFmt <- 'countdxw1'
    } else if ((SampleID == 'NESARCw2') | (SampleID == 'NCSR')) { # wave 2!
      varNameFmt <- 'countdxw2'
    }
    # Compute the variables that need to be added up
    vars <- ModelSpecs[[SampleID]][[ModelID]][['mmVars']]
    vars <- vars[!(vars %in% KeyVars[[SampleID]])]
    # Now, get me that overall variable!
    # Create matrix to temporarily store stuff
    m <- matrix(NA
                ,nrow=nrow(d)
                ,ncol=1)
    # Compute my stuff!
    # First, harmonize the variable names...
    colnames(m) <- varNameFmt
    # Now, get us the H indices!
    m[,1] <- d %>% select(all_of(vars)) %>% apply(1,sum)
    facCountVars <- as.data.frame(m)
    rv <- cbind(d,facCountVars)
  }
  return(rv)
}

# Compute factor split variables ------------------------------------------

# This function is to make the Script generator's life easier.
# Request type:
#   - 'overall' : overall dx count variable name
#   - 'facs' : factor dx count variable names
#   - 'dx' : variable names for all dx in the model
FacSplitNames <- function(SampleID,ModelID,ModelSpecs
                          ,KeyVars,RequestType) {
  rv <- NULL
  if (RequestType == 'overall') {
    if ((SampleID == 'NESARCw1') | (SampleID == 'NCS')) {
      rv <- 'countdxw1'
    } else if ((SampleID == 'NESARCw2') | (SampleID == 'NCSR')) {
      rv <- 'countdxw2'
    } 
  } else if ((RequestType == 'facs') & (ModelID != '1F')) {
    values <- names(MODEL_SPECS[[SampleID]][[ModelID]][["facs"]])
    # If bifactor model, remove G
    if(str_detect(ModelID,"Bi.")) {
      values <- values[values != 'G']
    }
    if ((SampleID == 'NESARCw1') | (SampleID == 'NCS')) {
      rv <- sprintf('count%sw1',values)
    } else if ((SampleID == 'NESARCw2') | (SampleID == 'NCSR')) {
      rv <- sprintf('count%sw2',values)
    } 
  } else if (RequestType == 'dx') {
    allvars <- MODEL_SPECS[[SampleID]][[ModelID]][['mmVars']]
    rv <- allvars[!(allvars %in% KeyVars[[SampleID]])]
  }
  return(rv)
}

computeFacSplit <- function(d # dataset
                            ,SampleID # sample ID
                            ,ModelID # model being run
                            ,ModelSpecs # from global
                            ,Crosswalks # from global
                            ,valid.Inputs # from global
) {
  rv <- NULL
  valid.Samples <- names(valid.Inputs)
  if (!(SampleID %in% valid.Samples)) {
    print(sprintf('INVALID SAMPLE SPECIFIED: %s! QUITTING.',SampleID))
  } else if ((!(ModelID %in% valid.Inputs[[SampleID]])) | (ModelID == '1F')) {
    print(sprintf('INVALID MODEL SPECIFIED FOR SAMPLE %s: %s',SampleID,ModelID))
    print(sprintf('VALID MODELS FOR THIS SAMPLE ARE %s'
                  ,paste(valid.Inputs[[SampleID]],collapse=' ')))
    print('QUITTING.')
  } else {
    # Variable names are formulaic.
    if ((SampleID == 'NESARCw1') | (SampleID == 'NCS')) {
      varNameFmt <- 'count%sw1'
    } else if ((SampleID == 'NESARCw2') | (SampleID == 'NCSR')) { # wave 2!
      varNameFmt <- 'count%sw2'
    }
    # Compute the variables that need to be added up
    # from the crosswalks
    facs <- ModelSpecs[[SampleID]][[ModelID]][['facs']]
    cw <- Crosswalks[[SampleID]]
    svars <- list()
    for(fac in names(facs)) {
      svars[[fac]] <- unname(
        as_vector(
          cw[cw$mylabs %in% facs[[fac]],'oglabs']
        )
      )
    }
    # Now, get me those splitting variables!
    # Create matrix to temporarily store stuff
    m <- matrix(NA
                ,nrow=nrow(d)
                ,ncol=length(svars))
    # Compute my stuff!
    # First, harmonize the variable names...
    colnames(m) <- sprintf(varNameFmt,names(svars))
    names(svars) <- sprintf(varNameFmt,names(svars))
    # Now, get us the H indices!
    for (facCountVar in names(svars)) {
      m[,facCountVar] <- d %>% select(svars[[facCountVar]]) %>% apply(1,sum)
    }
    facCountVars <- as.data.frame(m)
    rv <- cbind(d,facCountVars)
  }
  return(rv)
}

# createContainerDF -------------------------------------------------------

# This section will have the necessary functions to create a dataframe
# where Mplus outputs will be scraped. These functions are meant for the
# first stages, i.e., when only the measurement model is being examined.

# This function does most of the legwork, and works for a variety of
# model types. All that is necessary is:
#   - Simulations dataframne
#   - facs: named list of factors and corresponding diagnoses 
#           that load on them
            # Create a named list where:
            # Name = ONE capital letter for that factor
            # Elements = vector of all manifest variables to include
#   - phis: vector of correlations to be scraped for a given model.
            # Indicate the freely estimated factor correlations 
            # (phi) to be scraped (as the conjunction of capital letters)
createContainerDF_Helper <- function(Simulations,facs,phis,TableType,doDifftest,
                                     scrape_fac_variance) {
  # Get me those sim_IDs!
  sids <- Simulations$sim_ID
  # One row for each simulation
  n_rows <- length(sids)
  if (TableType == "EFA") {
    num_vars <- unlist(facs) %>% unique() %>% length()
    colnames <- c("errs.EFA","warn.EFA"
                  ,sprintf("Eigen_%02d",1:num_vars)
    )
    # Construct data frame
    n_cols <- length(colnames)
  } else {
    # Define some helpful variables
    fac_names <- names(facs) # list of abbreviated factor names
    n_facs <- length(fac_names) # Number of factors
    seprtr = '_' # Separator to use
    
    # Errors & model fit columns
    errmf <- c('errs','warn','N_obs','m_cfi','m_tli'
               ,'m_rmsea','m_srmr'
               ,'m_chi2','m_chi2df','m_chi2pv'
               ,'m_rmsea_90ci_lb','m_rmsea_90ci_ub'
               ,'m_rmsea_plt05'
    )
    
    if (doDifftest) {
      errmf <- c(errmf,'d_chi2','d_chi2df','d_chi2pv')
    }
    
    if (scrape_fac_variance) {
      fs <- names(facs)
      factor_variance_vars <- c(sprintf("var_%s",fs), # factor variance
                                sprintf("varse_%s",fs), # factor variance SE
                                sprintf("varpv_%s",fs)) # factor variance p-val
    }
    
    # Prefixes to combine for each factor
    prefixes <- c('l%s' # loading on factor (UNSTANDARDIZED)
                  ,'t%s' # loading on factor (STANDARDIZED)
                  ,'se%s' # SE of loading on factor (STANDARDIZED)
                  ,'pv%s' # p-val of loading on factor (STANDARDIZED)
    )
    
    # Compile all the relevant loading variables
    mmv <- c()
    vars <- c()
    for (lat.factor in fac_names) {
      tmp.prefixes <- sapply(prefixes,sprintf,lat.factor) 
      tmp.vars <- sapply(tmp.prefixes,paste,facs[[lat.factor]],sep=seprtr)
      vars <- c(vars,tmp.vars)
      mmv <- c(mmv,facs[[lat.factor]])
    }
    
    # Add now the stuff for the unstandardized residuals
    # get all the mm vars... manually because i'm lazy
    mmv <- unique(mmv)
    rv.prefixes <- c('re_%s' # loading on factor (UNSTANDARDIZED)
                  ,'rse_%s' # SE of loading on factor (STANDARDIZED)
                  ,'rpv_%s' # p-val of loading on factor (STANDARDIZED)
    )
    tmp.resvars <- as.vector(sapply(rv.prefixes,sprintf,mmv))
    vars <- c(vars,tmp.resvars)
    
    # R-square (r2cfa)
    r2.prefixes <- c('r2ce_%s' # rsquare (r2)
                     ,'r2cr_%s' # residual variance (1-r2)
                     )
    tmp.r2vars <- as.vector(sapply(r2.prefixes,sprintf,mmv))
    vars <- c(vars,tmp.r2vars)

    # Correlation data to scrape
    phi.prefixes <- c('phiu_%s' # phi correlation between two factors (UNSTANDARDIZED)
                      ,'phis_%s' # phi correlation between two factors (STANDARDIZED)
                      ,'phie_%s' # SE of phi correlation between two factors (STANDARDIZED)
                      ,'phip_%s' # p-val of phi correlation between two factors (STANDARDIZED)
    )
    
    # Add correlation data to list of variables if specified
    if ((!(is.na(phis))) && (length(phis) > 0)) {
      tmp.phivars <- as.vector(sapply(phi.prefixes,sprintf,phis))
      vars <- c(vars,tmp.phivars)
    }
    
    # add factor variance vars if specified
    if (scrape_fac_variance) {
      vars <- c(vars,factor_variance_vars)
    }
     
    colnames <- c(errmf,vars)
    
    # Construct data frame
    n_cols <- length(colnames)
  }
  # This is common to all.
  matx <- matrix(data = NA
                 ,nrow = n_rows
                 ,ncol = n_cols
                 ,dimnames = list(NULL
                                  ,colnames)
  )
  
  # Compile all these, PLUS Sim IDs, into a neat data frame!
  SimResults <- data.frame(sids,matx
                           ,stringsAsFactors = FALSE
  )
  # aaaand that's the output of our function. A blank dataframe!
  return(SimResults)
}

# createContainerDF: function that will generate an empty
# datafranme for storing scraping outputs. Note that this
# function will point you to the appropraite function for your
# model and sample.
#
#     Inputs: Simulations df
#             sample ID (NESARCw1, NESARCw2, NCS, NCSR)
#             model ID (2F, 3F, Bi2, Bi3)
#
# ADDITIONAL NOTE: For now, NESARC Wave 2 containers will look 
# identical to those from Wave 1. We will thus create wrapper 
# functions to use that point to each of these.
#
# TODO: update this function when I have the NCS functions.
#
createContainerDF <- function(Simulations,Sample,Model
                              ,valid.Inputs,ModelSpecs,TableType = "CFA"
                              ,doDifftest = FALSE) {
  rv <- NULL
  valid.Samples <- names(valid.Inputs)
  if (!(Sample %in% valid.Samples)) {
    print(sprintf('INVALID SAMPLE SPECIFIED: %s! QUITTING.',Sample))
  } else if (!(Model %in% valid.Inputs[[Sample]])) {
    print(sprintf('INVALID MODEL SPECIFIED FOR SAMPLE %s: %s',Sample,Model))
    print(sprintf('VALID MODELS FOR THIS SAMPLE ARE %s'
                  ,paste(valid.Inputs[[Sample]],collapse=' ')))
    print('QUITTING.')
  } else {
    # It should be the same structure for NESARC waves 1 or 2.
    facs <- ModelSpecs[[Sample]][[Model]][['facs']]
    phis <- ModelSpecs[[Sample]][[Model]]$phis
    
    scrape_fac_variance <- FALSE
    if (str_detect(Model,"vf")) {
      scrape_fac_variance <- TRUE
    }
    
    rv <- createContainerDF_Helper(Simulations,facs,phis,TableType,doDifftest,
                                   scrape_fac_variance)
  }
  return(rv)
}

# scrapeSim ---------------------------------------------------------------

# This section contains code for scraping simulation outputs for the
# measurement model into the dataframes created by createContainerDF.

getTargetLoad <- function(st_row,crosswalk,isSTDXY) {
  # NOTE: isSTDXY: TRUE if standardized, FALSE for unstandardized
  # This is assuming certain things (see functions above)
  # Factor
  fac <- str_sub(st_row['paramHeader'],start=1,end=1)
  # parameter
  pmt <- ifelse(st_row['parameter'] == 'est'
                ,ifelse(isSTDXY,'t','l')
                ,ifelse(isSTDXY
                        ,str_sub(st_row['parameter'],start=1,end=2)
                        ,NA)
  )
  # Dx
  dx <- crosswalk[crosswalk$oglabs == st_row['param'],'mylabs']
  # Label target
  tgt <- ifelse(is.na(pmt),NA,sprintf('%s%s_%s',pmt,fac,dx))
  return(tgt)
}

getTargetThreshold <- function(st_row,crosswalk) {
  # NOTE: isSTDXY: TRUE if standardized, FALSE for unstandardized
  # This is assuming certain things (see functions above)
  # Factor
  fac <- str_sub(st_row['paramHeader'],start=1,end=1)
  # parameter
  pmt <- ifelse(st_row['parameter'] == 'est'
                ,'e'
                ,str_sub(st_row['parameter'],start=1,end=2)
                )
  # Dx
  dx <- crosswalk[crosswalk$oglabs == st_row['param'],'mylabs']
  # Label target
  tgt <- ifelse(is.na(pmt),NA,sprintf('r%s_%s',pmt,dx))
  return(tgt)
}

getTargetR2 <- function(param,parameter,crosswalk) {
  # NOTE: isSTDXY: TRUE if standardized, FALSE for unstandardized
  # This is assuming certain things (see functions above)
  # parameter
  pmt <- ifelse(parameter == 'est'
                ,'e'
                ,'r')
  # Dx
  dx <- crosswalk$mylabs[which(crosswalk$oglabs == param)]
  # Label target
  tgt <- ifelse(is.na(pmt),NA,sprintf('r2c%s_%s',pmt,dx))
  return(tgt)
}


getTargetPhi <- function(st_row,phis,isSTDXY) {
  # NOTE: isSTDXY: TRUE if standardized, FALSE for unstandardized
  # This is assuming certain things (see functions above)
  
  # Factor labels
  fac1 <- str_sub(st_row['paramHeader'],start=1,end=1)
  fac2 <- str_sub(st_row['param'],start=1,end=1)
  # Identify the appropriate label based on model specification
  # E.g., if fac1="I" and fac2="E", label could be "IE" or "EI"
  # Want this to be consistent!
  rxp <- sprintf("((%s%s)|(%s%s))",fac1,fac2,fac2,fac1)
  label <- str_subset(phis,rxp) 
  
  # parameter
  pmt <- ifelse(st_row['parameter'] == 'est'
                ,ifelse(isSTDXY
                        ,'phis'
                        ,'phiu')
                ,ifelse(isSTDXY
                        ,ifelse(st_row['parameter'] == 'se'
                                ,'phie'
                                ,'phip')
                        ,NA)
  )
  # Label target
  tgt <- ifelse(is.na(pmt),NA,sprintf('%s_%s',pmt,label))
  return(tgt)
}

getTargetFacVariance <- function(factor, parameter) {
  if (parameter == 'est') {
    rv <- sprintf("var_%s",factor)
  } else if (parameter == 'se') {
    rv <- sprintf("varse_%s",factor)
  } else if (parameter == 'pval') {
    rv <- sprintf("varpv_%s",factor)
  } else { # should fix this later and use try catch, but UGHHHHH not in the mood
    rv <- NULL
  }
  rv
}

# This function will take care of scraping outputs from a CFA.
scrapeSim_Helper <- function(SimResults,SimPrefix,crosswalk_dx,facs,phis,doDifftest
                             ,scrape_fac_variance) {
  # Compute expected output filename
  OutFilename <- sprintf('Sim_%s-%s.out',SimPrefix,SimResults[1,'sids'])
  # Skip outputs that were already analyzed.
  if (!(is.na(SimResults[1,'errs']))) {
    print(sprintf('WARNING: Output %s was already scraped. Skipping.'
                  ,OutFilename))
  } else if (!(file.exists(OutFilename))){
    # The file does not exist. Skip over and print an error.
    print(sprintf('ERROR: Expected output file %s does not exist. Skipping.'
                  ,OutFilename))
  } else {
    # Read Mplus outputs
    sim_model <- readModels(OutFilename)
    # If an Mplus error was found, record it and leave the 
    # rest of the row values as NAs.
    if (!(length(sim_model$errors) == 0)) {
      errmsg <- unlist(sim_model$errors) %>% paste(collapse=" // ")
    } else {
      errmsg <- ""
    }
    
    if (errmsg != "" & !(startsWith(errmsg,"FACTOR SCORES COULD NOT BE COMPUTED."))) {
      SimResults[1,'errs'] <- errmsg
    } else {
      # If we've made it here, then the model was estimated. Yay!
      
      # save factor scores if possible
      if("savedata" %in% names(sim_model)) {
        fs_filename <- sprintf('Sim_%s-%s-FS.csv',SimPrefix,SimResults[1,'sids'])
        write_csv(sim_model$savedata,fs_filename,na=".")
      }
      
      # Let's start off by extracting measurement model & regression parameters
      sim_StdXY_Load_params <- paramExtract(sim_model$parameters$stdyx.standardized,
                                            params = 'loading')
      sim_UnStd_Load_params <- paramExtract(sim_model$parameters$unstandardized,
                                            params = 'loading')
      sim_UnStd_Thr <- paramExtract(sim_model$parameters$unstandardized,
                                    params = 'e')
      
      
      cdx <- quo(crosswalk_dx)
      
      sim_r2 <- sim_model$parameters$r2 %>%
                  select(param,est,resid_var) %>%
                  gather(est,resid_var
                         ,key="parameter",value="value") %>%
                  mutate(target = map2_chr(param,parameter
                                           , ~ getTargetR2(.x,.y
                                                           ,!!cdx))) %>%
                  select(target,value)
      
      if (scrape_fac_variance) {
        sim_fvar <- sim_model$parameters$unstandardized %>%
          paramExtract(params='variability') %>%
          select(param,est,se,pval) %>%
          gather(est,se,pval
                 ,key="parameter",value="value") %>%
          rename(factor=param) %>%
          mutate(target = map2_chr(factor,parameter
                                   ,getTargetFacVariance)) %>%
          select(target,value)
      }
      
      # For thresholds, we only want thresholds and not means
      sim_UnStd_Thr <- sim_UnStd_Thr[sim_UnStd_Thr$paramHeader == "Thresholds",]
      # Clean up param and reove that weird $1
      sim_UnStd_Thr$param <- gsub("\\$1","",sim_UnStd_Thr$param)
      
      # We don't care about est_se
      sLp_StdXY_abbr <- select(sim_StdXY_Load_params,-est_se)
      sLp_UnStd_abbr <- select(sim_UnStd_Load_params,-est_se)
      sTp_UnStd_abbr <- select(sim_UnStd_Thr,-est_se)

      # Convert to wide format
      sLp_StdXY_long <- gather(sLp_StdXY_abbr
                               ,est,se,pval
                               ,key="parameter"
                               ,value='value')
      sLp_UnStd_long <- gather(sLp_UnStd_abbr
                               ,est,se,pval
                               ,key="parameter"
                               ,value='value')
      sTp_UnStd_long <- gather(sTp_UnStd_abbr
                               ,est,se,pval
                               ,key="parameter"
                               ,value='value')
      # Compute target
      sLp_StdXY_long$target <- apply(sLp_StdXY_long,1,getTargetLoad
                                     ,crosswalk_dx, TRUE) 
      sLp_UnStd_long$target <- apply(sLp_UnStd_long,1,getTargetLoad
                                     ,crosswalk_dx, FALSE) 
      sTp_UnStd_long$target <- apply(sTp_UnStd_long,1,getTargetThreshold
                                     ,crosswalk_dx) 
      
      sLp_StdXY_longAbr <- select(sLp_StdXY_long,target,value)
      sLp_UnStd_longAbr <- select(sLp_UnStd_long,target,value)
      sLp_UnStd_longAbr <- sLp_UnStd_longAbr[!(is.na(sLp_UnStd_longAbr$target)),]
      
      sTp_UnStd_longAbr <- select(sTp_UnStd_long,target,value)
      
      sLp_longAbr <- rbind(sLp_StdXY_longAbr,sLp_UnStd_longAbr
                           ,sTp_UnStd_longAbr,sim_r2)
      if (scrape_fac_variance) {
        sLp_longAbr <- rbind(sLp_longAbr,sim_fvar)
      }
      
      # This next part does the same, but for pertinent correlations
      # (if applicable!)
      if (!(is.na(phis))) {
        # Get factor correlation data
        sim_UndU_params <- paramExtract(sim_model$parameters$unstandardized,
                                       params = 'undirected')
        sim_UndS_params <- paramExtract(sim_model$parameters$stdyx.standardized,
                                        params = 'undirected')
        # We don't care about est_se
        sim_UndU_params_abbr <- select(sim_UndU_params,-est_se)
        sim_UndS_params_abbr <- select(sim_UndS_params,-est_se)
        # Convert to wide format
        sim_UndU_params_long <- gather(sim_UndU_params_abbr
                                       ,est,se,pval
                                       ,key="parameter"
                                       ,value='value')
        sim_UndS_params_long <- gather(sim_UndS_params_abbr
                                       ,est,se,pval
                                       ,key="parameter"
                                       ,value='value')
        # Compute target
        sim_UndS_params_long$target <- apply(sim_UndS_params_long,1,getTargetPhi
                                             ,phis, TRUE) 
        sim_UndU_params_long$target <- apply(sim_UndU_params_long,1,getTargetPhi
                                             ,phis, FALSE) 
        
        UndS_longAbr <- select(sim_UndS_params_long,target,value)
        UndU_longAbr <- select(sim_UndU_params_long,target,value)
        UndU_longAbr <- UndU_longAbr[!(is.na(UndU_longAbr$target)),]
        
        sLp_longAbr <- rbind(sLp_longAbr,UndS_longAbr,UndU_longAbr)
        
      }
      
      # Throw in the stuff we want
      # This is hack-ish, but should work...
      tr <- spread(sLp_longAbr,key=c('target'),value=c('value')) %>% as_tibble()
      row <- SimResults[1,] %>% as_tibble() %>% bind_rows(tr)
      row$sids[2] <- row$sids[1]
      row <- row[2,]
      SimResults[1,] <- row
      if (length(sim_model$errors) != 0) {
        SimResults[1,'errs'] <- errmsg
      }
      # Model fit parameters
      SimResults[1,'m_cfi'] <- sim_model$summaries$CFI
      SimResults[1,'m_tli'] <- sim_model$summaries$TLI
      SimResults[1,'m_rmsea'] <- sim_model$summaries$RMSEA_Estimate
      SimResults[1,'m_rmsea_90ci_lb'] <- sim_model$summaries$RMSEA_90CI_LB
      SimResults[1,'m_rmsea_90ci_ub'] <- sim_model$summaries$RMSEA_90CI_UB
      SimResults[1,'m_srmr'] <- sim_model$summaries$SRMR
      SimResults[1,'m_chi2'] <- sim_model$summaries$ChiSqM_Value
      SimResults[1,'m_chi2df'] <- sim_model$summaries$ChiSqM_DF
      SimResults[1,'m_chi2pv'] <- sim_model$summaries$ChiSqM_PValue
      # N observations
      SimResults[1,'N_obs'] <- sim_model$summaries$Observations
      # Warnings
      if (length(sim_model$warnings) > 0) {
        SimResults[1,'warn'] <- unlist(sim_model$warnings) %>% paste(collapse=" // ")
      }
      if (doDifftest) {
        if ("ChiSqDiffTest_Value" %in% names(sim_model$summaries)) {
          SimResults[1,'d_chi2'] <- sim_model$summaries$ChiSqDiffTest_Value
        }
        if ("ChiSqDiffTest_DF" %in% names(sim_model$summaries)) {
          SimResults[1,'d_chi2df'] <- sim_model$summaries$ChiSqDiffTest_DF
        }
        if ("ChiSqDiffTest_PValue" %in% names(sim_model$summaries)) {
          SimResults[1,'d_chi2pv'] <- sim_model$summaries$ChiSqDiffTest_PValue
        }
      }
    }
  }
  
  return(SimResults)
}

scrapeSim <- function( SimResults # this is actually just a row, not the full thing!
                      ,SimPrefix
                      ,Sample
                      ,Model
                      ,valid.Inputs
                      ,Crosswalks
                      ,ModelSpecs
                      ,doDifftest = FALSE
                      ) {
    valid.Samples <- names(valid.Inputs)
    if (!(Sample %in% valid.Samples)) {
      print(sprintf('INVALID SAMPLE SPECIFIED: %s! QUITTING.',Sample))
    } else if (!(Model %in% valid.Inputs[[Sample]])) {
      print(sprintf('INVALID MODEL SPECIFIED FOR SAMPLE %s: %s',Sample,Model))
      print(sprintf('VALID MODELS FOR THIS SAMPLE ARE %s'
                    ,paste(valid.Inputs[[Sample]],collapse=' ')))
      print('QUITTING.')
    } else {
      facs <- ModelSpecs[[Sample]][[Model]][['facs']]
      phis <- ModelSpecs[[Sample]][[Model]][['phis']]
      mmVars <- ModelSpecs[[Sample]][[Model]][['mmVars']]
      crosswalk <- Crosswalks[[Sample]]
      crosswalk <- crosswalk[(crosswalk$oglabs %in% mmVars),]
      scrape_fac_variance <- ifelse(str_detect(Model,"vf"),TRUE,FALSE)
      SimResults <- scrapeSim_Helper(SimResults
                                     ,SimPrefix
                                     ,crosswalk
                                     ,facs
                                     ,phis
                                     ,doDifftest
                                     ,scrape_fac_variance)
    }
    return(SimResults)
}

scrapeSim.EFA <- function (SimResults,SimPrefix) {
  # Compute expected output filename
  OutFilename <- sprintf('Sim_%s-%s_EFA.out',SimPrefix,SimResults[1,'sids'])
  if (!(file.exists(OutFilename))){
    # The file does not exist. Skip over and print an error.
    print(sprintf('ERROR: Expected output file %s does not exist. Skipping.'
                  ,OutFilename))
  } else {
    # Read Mplus outputs
    sim_model <- readModels(OutFilename)
  
    # If an Mplus error was found, record it and leave the 
    # rest of the row values as NAs.
    if (!(length(sim_model$errors) == 0)) {
      # Save our loooong error list
      SimResults$errs.EFA[1] <- unlist(sim_model$errors) %>% paste(collapse = " // ")
    } else {
      SimResults$errs.EFA[1] <- NA
    }
    
    if (!(length(sim_model$warnings) == 0)) {
      # Save our loooong error list
      SimResults$warn.EFA[1] <- unlist(sim_model$warnings) %>% paste(collapse = " // ")
    } else {
      SimResults$warn.EFA[1] <- NA
    }
    
    # Check if exists
    if ("gh5" %in% names(sim_model)) {
      efa_eigen <- as.numeric(sim_model$gh5$efa$eigenvalues)
      # Check how many Eigen columns we created previously.
      n_eigen <- sum(str_detect(names(SimResults),"Eigen_"))
      # Coerce to get the number of Eigenvalues expected.
      if (length(efa_eigen) < n_eigen) {
        efa_eigen <- c(efa_eigen,rep(NA,n_eigen-length(efa_eigen)))
      } else if (length(efa_eigen) > n_eigen) {
        efa_eigen <- efa_eigen[1:n_eigen]
      }
      efa_eigen.m <- matrix(efa_eigen,ncol=n_eigen,nrow=1)
      colnames(efa_eigen.m) <- sprintf("Eigen_%02d",1:n_eigen)
      tr <- as_tibble(efa_eigen.m)
      SimResults[1,names(tr)] <- tr
    }
  }
  return(SimResults)
}

# H-Index Calculation -----------------------------------------------------

# For a given vector of Lambdas (loadings)
H.index.formula <- function(lambdas_v) {
  lambdas_v_sq <- lambdas_v^2
  H <- (1/(
    1 + (1 /
           sum(lambdas_v_sq / (1 - lambdas_v_sq))
    ) )
  )
  return(H)
}

getFS <- function(ind,Sims,l) {
  sim_result <- as_tibble(cbind(name=names(Sims),t(Sims[ind,]))) %>%
    rename_if(!startsWith(names(.),"name"),~"value")
  
  sr_q <- quo(sim_result)
  
  fs <- enframe(l) %>%
    mutate(
      loadings = map(value, ~ !!sr_q %>% filter(name %in% .x) 
                     %>% pull(value) 
                     %>% unlist()
                     %>% as.numeric()
                     )
      ,valid = map_lgl(loadings,~ifelse((length(.x) - sum(is.na(.x))) != length(.x)
                                          ,FALSE,TRUE))
      ,sum_load = map2_dbl(loadings,valid, ~ifelse(.y, sum(.x), NA))
      ,sum_of_sq_load = map2_dbl(loadings,valid, ~ifelse(.y, sum(.x^2), NA))
      ,sum_load_sq = sum_load^2
    )
  
  return(fs)
}

compute_omega_row <- function(ind,Sims,l) {
  fs <- getFS(ind,Sims,l)
  
  if (fs$valid[which(fs$name == "ResidVar")]) {
    numerator <- fs %>% 
                  filter(name != "ResidVar") %>% 
                  pull(sum_load_sq) %>%
                  sum()
    denominator <- fs %>% 
                    filter(name == "ResidVar") %>%
                    pull(sum_load) %>% 
                    sum(., numerator)
    
    stopifnot(length(numerator) == 1,length(denominator) == 1)
    
    omega <- numerator / denominator
  } else {
    omega <- NA
  }
  
  return(omega)
}

compute_omegaH_row <- function(ind,Sims,l,factr) {
  fs <- getFS(ind,Sims,l)
  stopifnot(factr %in% fs$name)
  
  if (fs$valid[which(fs$name == "ResidVar")]) {
    
    numerator <- fs %>% 
      filter(name == factr) %>% 
      pull(sum_load_sq) 
    
    denom_rv <- fs %>% 
      filter(name == "ResidVar") %>%
      pull(sum_load)
    
    denominator <- fs %>%
      filter(name != "ResidVar") %>%
      pull(sum_load_sq) %>%
      sum(.,denom_rv)
    
    stopifnot(length(numerator) == 1,length(denominator) == 1)
    
    omegaH <- numerator / denominator
  } else {
    omegaH <- NA
  }
  
  return(omegaH)
}


compute_ecv_row <- function(ind,Sims,l) {
  fs <- getFS(ind,Sims,l)
  stopifnot("G" %in% fs$name)
  
  if (fs$valid[which(fs$name == "ResidVar")]) {
    numerator <- fs %>%
      filter(name == "G") %>%
      pull(sum_of_sq_load)
    
    denominator <- fs %>%
      filter(name != "ResidVar") %>%
      pull(sum_of_sq_load) %>%
      sum()
    
    stopifnot(length(numerator) == 1,length(denominator) == 1)
    
    ecv <- numerator / denominator
  } else {
    ecv <- NA
  }
  
  return(ecv)
}

compute_omegaHS_row <- function(ind,Sims,l,subscale) {
  fs <- getFS(ind,Sims,l)
  stopifnot("G" %in% fs$name)
  stopifnot(subscale %in% fs$name)
  
  if (fs$valid[which(fs$name == "ResidVar")]) {
  
    numerator <- fs %>% 
      filter(name == subscale) %>% 
      pull(sum_load_sq) 
    
    denom_rv <- fs %>% 
      filter(name == "ResidVar") %>%
      pull(sum_load)
    
    denominator <- fs %>%
      filter(name %in% c("G",subscale)) %>%
      pull(sum_load_sq) %>%
      sum(.,denom_rv)
    
    stopifnot(length(numerator) == 1,length(denominator) == 1)
    
    omegaHS <- numerator / denominator
  } else {
    omegaHS <- NA
  }
  
  return(omegaHS)
}

compute_Hindex <- function(Sims,facs) {
  fac_names <- names(facs)
  prefixes <- sprintf("t%s_",fac_names)
  
  # Determine target variable names for each factor (STANDARDIZED loadings)
  l <- mapply(function(dxs, facLabel) {sprintf("t%s_%s",facLabel,dxs)}
              ,facs,names(facs),SIMPLIFY = FALSE)
  
  ind <- 1:nrow(Sims)
  
  # preallocate matrix
  targets <- sprintf("Hindx_%s",fac_names)
  mat <- matrix(data=NA
                ,nrow=nrow(Sims)
                ,ncol=length(targets)
                ,dimnames = list(NULL,targets))
  
  for (factor in fac_names) {
    target <- sprintf("Hindx_%s",factor)
    for (i in ind) {
      lambdas <- Sims[i,] %>% 
                  select(l[[factor]]) %>% 
                  as.numeric()
      if (length(lambdas) == (length(lambdas) - sum(is.na(lambdas)))) {
        hind <- H.index.formula(lambdas)
        mat[i,target] <- hind
      }
    }
  }
  
  rv <- bind_cols(Sims,as_tibble(mat))
  
  return(rv)
}

compute_omega <- function(Sims,facs) {
  fac_names <- names(facs)
  
  # Determine target variable names for each factor (STANDARDIZED loadings)
  l <- mapply(function(dxs, facLabel) {sprintf("t%s_%s",facLabel,dxs)}
              ,facs,names(facs),SIMPLIFY = FALSE)
  # add residuals
  l$ResidVar <- unlist(facs) %>% 
                  unique() %>% 
                  map_chr(~sprintf("r2cr_%s",.x))
  
  ind <- 1:nrow(Sims)
  rv <- Sims
  omega <- map_dbl(ind, ~ compute_omega_row(.x,Sims,l))
  return(omega)
}

compute_ecv <- function(Sims,facs) {
  fac_names <- names(facs)
  
  # Determine target variable names for each factor (STANDARDIZED loadings)
  l <- mapply(function(dxs, facLabel) {sprintf("t%s_%s",facLabel,dxs)}
              ,facs,names(facs),SIMPLIFY = FALSE)
  # add residuals
  l$ResidVar <- unlist(facs) %>% 
    unique() %>% 
    map_chr(~sprintf("r2cr_%s",.x))
  
  ind <- 1:nrow(Sims)
  rv <- Sims
  ecv <- map_dbl(ind, ~ compute_ecv_row(.x,Sims,l))
  return(ecv)
}

# compute_omegaH <- function(Sims,facs) {
#   fac_names <- names(facs)
#   prefixes <- sprintf("")
#   
#   # Determine target variable names for each factor (STANDARDIZED loadings)
#   l <- mapply(function(dxs, facLabel) {sprintf("t%s_%s",facLabel,dxs)}
#               ,facs,names(facs),SIMPLIFY = FALSE)
#   # add residuals
#   l$ResidVar <- unlist(facs) %>% 
#     unique() %>% 
#     map_chr(~sprintf("r2cr_%s",.x))
#   
#   ind <- 1:nrow(Sims)
#   rv <- Sims
#   omegaH <- map_dbl(ind, ~ compute_omegaH_row(.x,Sims,l))
#   return(omegaH)
# }

compute_omegaH_relomega <- function(Sims,facs) {
  fac_names <- names(facs)
  
  # Determine target variable names for each factor (STANDARDIZED loadings)
  l <- mapply(function(dxs, facLabel) {sprintf("t%s_%s",facLabel,dxs)}
              ,facs,names(facs),SIMPLIFY = FALSE)
  # add residuals
  l$ResidVar <- unlist(facs) %>% 
    unique() %>% 
    map_chr(~sprintf("r2cr_%s",.x))
  ind <- 1:nrow(Sims)
  omegaH_relomega_list <- list()
  for (factr in fac_names) {
    target_omegaH <- sprintf("omegaH_%s",factr)
    target_omegaRel <- sprintf("relOmega_%s",factr)
    tmp_omegaH <- map_dbl(ind, ~ compute_omegaH_row(.x,Sims,l,factr))
    tmp_omegaRel <- tmp_omegaH / Sims$omega
    omegaH_relomega_list[[target_omegaH]] <- tmp_omegaH
    omegaH_relomega_list[[target_omegaRel]] <- tmp_omegaRel
  }
  omegaH_rel_df <- as_tibble(omegaH_relomega_list)
  rv <- bind_cols(Sims,omegaH_rel_df)
  return(rv)
}


compute_omegaHS <- function(Sims,facs) {
  subscales <- names(facs)[names(facs) != "G"]
  fac_names <- names(facs)
  
  # Determine target variable names for each factor (STANDARDIZED loadings)
  l <- mapply(function(dxs, facLabel) {sprintf("t%s_%s",facLabel,dxs)}
              ,facs,names(facs),SIMPLIFY = FALSE)
  # add residuals
  l$ResidVar <- unlist(facs) %>% 
    unique() %>% 
    map_chr(~sprintf("r2cr_%s",.x))
  
  ind <- 1:nrow(Sims)
  omegaHS_list <- list()
  for (subscale in subscales) {
    target <- sprintf("omegaHS_%s",subscale)
    tmp.omegaHS <- map_dbl(ind, ~ compute_omegaHS_row(.x,Sims,l,subscale))
    omegaHS_list[[target]] <- tmp.omegaHS
  }
  omegaHS_df <- as_tibble(omegaHS_list)
  rv <- bind_cols(Sims,omegaHS_df)
  return(rv)
}

compute_puc <- function(facs) {
  # Basedon constantinou & Fonagy preprint formula
  indicators <- unlist(facs) %>% unique()
  p <- length(indicators)
  
  p_ms <- facs
  p_ms[['G']] <- NULL
  p_ms <- lapply(p_ms,length) %>% unlist() %>% as.numeric()
  numerator_elements <- p_ms*(p_ms-1)/2
  numerator <- sum(numerator_elements)
  denominator <- (p*(p-1)) / 2
  puc <- 1 - (numerator / denominator)
  return(puc)
}

# Minor helper functions --------------------------------------------------

# A little helper function to compile a character vector into
# a literal string representation of itself
cv_to_str <- function(cv) {
  rv <- str_c("c("
              ,str_c(sprintf("'%s'",cv)
                     ,collapse=",")
              ,")"
  )
  return(rv)
}



# This file has a bunch of info on the different samples we used, basically
# we used globals to reduce manual inputs across and have it in a centralized
# place...

# Implemented samples -----------------------------------------------------

VALID_SAMPLES <- c( "NESARCw1"
                   ,"NESARCw2"
                   ,"NCS"
                   ,"NCSR"
                   )

# Datasets ----------------------------------------------------------------

DATASET_NAMES <- list( "NESARCw1" = "NESARCw1w2_Processed_AllIDs.csv"
                       ,"NESARCw2" = "NESARCw1w2_Processed_W2only.csv"
                       ,"NCS" = "NCS_Processed_Pt2IDsOnly.csv"
                       ,"NCSR" = "NCSR_Processed_Pt2IDsOnly.csv"
)


# Key variables -----------------------------------------------------------

# Specify the ID, weight, and stratification variables that must ALWAYS
# be included in a model.

KEY_VARS <- list ('NESARCw1' = c(id='IDNUM',strat='W1STRAT',clust='W1PSU',wt='W1WT')
                  ,'NESARCw2' = c(id='IDNUM',strat='W2STRAT',clust='W2PSU',wt='W2WT')
                  ,'NCS' = c(id='CASEID',strat='STR',wt='P2TOBWT')
                  ,'NCSR' = c(id='CASEID',strat='SESTRAT',clust='SECLUSTR',wt='NCSRWTLG')
                  )

# Placeholder variables ---------------------------------------------------

CROSSWALKS <- list()
MODEL_SPECS <- list()


# External variable step preferences --------------------------------------

# Should regressions be computed using factor score estimates or by fixing
# the loadings of the measurement model and running a SEM model?
E.METHOD <- "BOTH" # Valid inputs: "FSE", "FIXED_SEM", "BOTH"

# Should the loadings of psychopathology factors on each individual 
# demographic variable be tested, simultaneously?
# E.g., test FACTOR on AGE SEX RACE
E.DO_FAC_DEMOS <- TRUE # Valid inputs: TRUE or FALSE

# Should we compute a set of analyses where demographics are added as covariates?
# E.g., if interested in HEALTH, test as FACTOR on HEALTH AGE SEX.
# IMPORANT: This MUST be marked true if you want to compute R-squared due to
#           psychopathology.
E.DO_E_FAC_DEMOS <- TRUE

# Should we compute a set of analyses where demographics are NOT included as
# covariates?
# E.g., if interested in HEALTH, test as FACTOR on HEALTH.
E.DO_E_FAC <- TRUE

# Should we compute a set of analyses where we regress each variable on the
# demographics only? 
# E.g., if interested in HEALTH, test as HEALTH on AGE SEX.
# IMPORANT: This MUST be marked true if you want to compute R-squared due to
#           psychopathology.
E.DO_E_DEMOS <- TRUE

# (FOR CORRELATED FACTOR MODELS ONLY)
# Include the other factor as a covariate in the models?
# E.g., in a 2F model (INT/EXT) where interested in HEALTH, test as
#       HEALTH on INT EXT
E.INCL_OTHER_FAC <- "BOTH" # Valid inputs: "INLCUDE", "EXCLUDE", "BOTH"
                           # "BOTH" will run it with and without the other fac

# Compile all these preferences into a neat little list (will be helpful when
# passing these guys on to a function)
EXT_PREFS <- list( "METHOD" = E.METHOD
                  ,"DO_FAC_DEMOS" = E.DO_FAC_DEMOS
                  ,"DO_E_FAC_DEMOS" = E.DO_E_FAC_DEMOS
                  ,"DO_E_FAC" = E.DO_E_FAC
                  ,"DO_E_DEMOS" = E.DO_E_DEMOS
                  ,"INCL_OTHER_FAC" = E.INCL_OTHER_FAC
)
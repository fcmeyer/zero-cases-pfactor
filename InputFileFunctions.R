# File has some helpful stuff for writing mplus input files... used some of 
# these two wrangle our various input files together.


# Helper functions --------------------------------------------------------

VarsToSplitLines <- function(vec.char
                             ,add.semicolon=TRUE # add semicolon on last line
                             ,add.star=FALSE
) {
  len <- length(vec.char)
  if (len < 6) {
    l <- list(vec.char)
    nlines <- 1
  } else {
    if (len %% 5 == 0) {
      nlines <- len / 5
      fac <- rep(1:nlines,each=5)
    } else {
      nlines <- floor(len/5) + 1
      fac <- c(rep(1:(nlines-1),each=5)
               ,rep(nlines,each=(len %%5)))
    }
    l <- split(vec.char,fac)
  }
  if(add.star) {
    l[[1]][1] <- str_c(l[[1]][1],"*",sep="")
  }
  l <- lapply(l,paste,sep="",collapse=" ")
  if (add.semicolon) {
    l[[nlines]] <- str_c(l[[nlines]],";",sep="")
  }
  return(l)
}


# Generic Prefix Functions ------------------------------------------------

L.CatVars <- function(catVars) {
  cv <- catVars
  L <- VarsToSplitLines(cv
                        ,add.semicolon=TRUE,add.star=FALSE)
  L <- c("CATEGORICAL ARE",L,"")
  return(L)
}

L.UseVars <- function(useVars) {
  uv <- useVars
  L <- VarsToSplitLines(uv
                        ,add.semicolon=TRUE,add.star=FALSE)
  L <- c("USEVARIABLES ARE",L,"")
  return(L)
}

L.IDSpecs <- function(sample,KEY_VARS) {
  idvars <- KEY_VARS[[sample]]
  
  id <- sprintf("IDVARIABLE IS %s;",KEY_VARS[[sample]][["id"]])
  strat <- sprintf("STRATIFICATION IS %s;",KEY_VARS[[sample]][["strat"]])
  wt <- sprintf("WEIGHT IS %s;",KEY_VARS[[sample]][["wt"]])
  
  rv <- NULL
  if ("clust" %in% names(idvars)) {
    clust <- sprintf("CLUSTER IS %s;",KEY_VARS[[sample]][["clust"]])
    rv <- list(id,strat,clust,wt)
  } else {
    rv <- list(id,strat,wt)
  }
  
  rv <- c(rv,"")
  
  return(rv)
}

L.Analysis <- function(estimator = "WLSMV") {
  rv <- list("ANALYSIS:"
             ,""
             ,"TYPE = COMPLEX;"
             ,sprintf("ESTIMATOR = %s;",estimator)
             ,""
  )
  return(rv)
}

# Generic Suffix Functions ------------------------------------------------

# If we want to save factor scores!
L.SaveFS <- function(SIM_PREFIX,sim_id) {
  rv <- list(""
             ,"SAVEDATA:"
             ,""
             ,sprintf("FILE = Sim_%s-%s.txt;",SIM_PREFIX,sim_id)
             ,"SAVE = FSCORES;"
             ,"")
  
  return(rv)
}

# CFA + EFA Wrappers ------------------------------------------------------

L.UseVars.MM <- function(sample,model
                         ,MODEL_SPECS,KEY_VARS) {
  vars <- MODEL_SPECS[[sample]][[model]][["mmVars"]]
  vars <- vars[!(vars %in% KEY_VARS[[sample]])]
  return(L.UseVars(vars))
}

L.CatVars.MM <- function(sample,model
                         ,MODEL_SPECS) {
  varTab <- MODEL_SPECS[[sample]][[model]][["varTypes"]]
  catVars <- varTab$oglab[varTab$type == "Categorical"]
  return(L.CatVars(catVars))
}

# CFA Functions -----------------------------------------------------------

L.ModelCFA <- function(sample,model
                       ,MODEL_SPECS,CROSSWALKS,KEY_VARS) {
  # Determine if we're dealing with a bifactor model
  # (note: all of them are done orthogonal)
  isBifactor <- str_detect(model,"Bi")
  
  # Create lines for the loadings
  
  # Identify variables used in measurement model
  mmVarsNID <- MODEL_SPECS[[sample]][[model]][["mmVars"]]
  mmVarsNID <- mmVarsNID[!(mmVarsNID %in% KEY_VARS[[sample]])]
  
  # Use crosswalk to determine names of vars loading on
  # each factor
  cw <- CROSSWALKS[[sample]]
  cw <- cw[cw$oglabs %in% mmVarsNID,]
  
  # Pull model specs
  facs <- MODEL_SPECS[[sample]][[model]][["facs"]]
  phis <- MODEL_SPECS[[sample]][[model]][["phis"]]
  
  # Create the FACTOR BY VAR* ... VAR; segment
  rv <- list("MODEL:","")
  for (factor in names(facs)) {
    vars <- cw$oglabs[cw$mylabs %in% facs[[factor]]]
    L <- VarsToSplitLines(vars,add.semicolon=TRUE
                          ,add.star=TRUE)
    L[[1]] <- str_c(sprintf("%s BY ",factor)
                    ,L[[1]])
    # Add blank line between the factors!
    L <- list.append(L,"")
    # Append that
    rv <- c(rv,L)
  }
  
  # Fix factor means and variances!
  for (factor in names(facs)) {
    line <- sprintf("[%s@0]; %s@1;",factor,factor)
    rv <- c(rv,line)
  }
  rv <- c(rv,"")
  
  # If we have more than one factor present, we should deal
  # with "WITH" statements.
  if (length(names(facs)) > 1) {
    # Multi factor models: set params free!
    fac_names <- names(facs)
    
    # Generate lines for correlations
    phi.mat <- combinations(fac_names,k=2,replace=FALSE)
    prev <- ""
    for (i in 1:nrow(phi.mat)) {
      line.tmp <- sprintf("%s WITH %s"
                          ,phi.mat[i,1],phi.mat[i,2])
      if (isBifactor) {
        out <- str_c(line.tmp,"@0;")
      } else if (prev == phi.mat[i,1]) {
        out <- str_c(line.tmp,";")
      } else {
        out <- str_c(line.tmp,"*;")
      }
      prev <- phi.mat[i,1]
      rv <- c(rv,out)
    }
  }
  rv <- c(rv,"")
  return(rv)
}

# Returns list of character vectors for write_lines()
AssembleFile.CFA <- function(sample,model
                             ,MODEL_SPECS,CROSSWALKS,KEY_VAR
                             ,difftest=FALSE) {
  rv <- list()
  
  # Add USEVARIABLES
  tmp <- L.UseVars.MM(sample,model
                      ,MODEL_SPECS,KEY_VARS)
  rv <- c(rv,tmp)
  
  # Add CATEGORICAL statement
  tmp <- L.CatVars.MM(sample,model
                      ,MODEL_SPECS)
  rv <- c(rv,tmp)
  
  # Add ID/STRAT/CLUSTER/WT section
  tmp <- L.IDSpecs(sample,KEY_VARS)
  rv <- c(rv,tmp)

  # Add ANALYSIS section
  tmp <- L.Analysis("WLSMV")
  rv <- c(rv,tmp)
  
  # Add difftest line if applicable
  if (difftest) {
    tmp <- c("DIFFTEST = GSUB_DIFFFILENAME;","")
    rv <- c(rv,tmp)
  }
  
  # Add MODEL section
  tmp <- L.ModelCFA(sample,model
                    ,MODEL_SPECS,CROSSWALKS,KEY_VARS)
  rv <- c(rv,tmp)
  
  # Add OUTPUT section
  tmp <- "OUTPUT: SAMPSTAT SVALUES CINT RESIDUAL TECH1 TECH4 STANDARDIZED;"
  rv <- c(rv,tmp,"")
  
  return(rv)
}


# EFA Functions -----------------------------------------------------------

# Includes analysis, plot, and output statement.
L.SuffixEFA <- function(numfacs) {
  rv <- list("ANALYSIS:"
             ,""
             ,"TYPE = COMPLEX;"
             ,sprintf("TYPE = EFA (1 %d);",numfacs)
             ,""
             ,"PLOT: TYPE = PLOT2;"
             ,""
             ,"OUTPUT: SAMPSTAT;"
  )
  return(rv)
}

# Returns list of character vectors for write_lines()
AssembleFile.EFA <- function(sample,model
                             ,MODEL_SPECS,CROSSWALKS,KEY_VAR) {
  rv <- list()
  
  # Add USEVARIABLES
  tmp <- L.UseVars.MM(sample,model
                      ,MODEL_SPECS,KEY_VARS)
  rv <- c(rv,tmp)
  
  # Add CATEGORICAL statement
  tmp <- L.CatVars.MM(sample,model
                      ,MODEL_SPECS)
  rv <- c(rv,tmp)
  
  # Add ID/STRAT/CLUSTER/WT section
  tmp <- L.IDSpecs(sample,KEY_VARS)
  rv <- c(rv,tmp)
  
  # Add ANALYSIS section
  tmp <- L.SuffixEFA(10)
  rv <- c(rv,tmp)
  
  return(rv)
}


# External Model Functions ------------------------------------------------

# Change the path to refer to the prior directory!
L.EditPrefixPath.EM <- function(MplusPrefix) {
  prefix <- MplusPrefix
  line.data <- prefix[[str_which(prefix,"DATA: ")]]
  prefix[[str_which(prefix,"DATA: ")]] <- str_replace(line.data
                                                      ,'FILE = \"'
                                                      ,'FILE = \"../')
  return(prefix)
}

L.ModelEM.RegrFull <- function(evmod) {
  y <- evmod$y
  xs <- evmod$xs
  line <- VarsToSplitLines(c(y,"ON",xs)
                           ,add.semicolon = TRUE
                           ,add.star = FALSE)
  rv <- list(line,"")
  return(rv)
}

L.ModelEM.FixLoad <- function(evmod,Sims,SimID) {
  rv <- list()
  for (factor in names(evmod$model_struct$facs)) {
    for (dx in evmod$model_struct$facs[[factor]]) {
      position <- match(SimID,Sims$sids)
      colname <- sprintf("l%s_%s",factor,dx)
      value <- pull(Sims,colname)[position]
      varname <- evmod$model_struct$varTypes$oglab[
        match(dx,evmod$model_struct$varTypes$mylab)
        ]
      tmp.line <- sprintf("%s BY %s@%f;",factor,varname,value)
      rv <- c(rv,tmp.line)
    }
    rv <- c(rv," ")
  }
  return(rv)
}

L.ModelEM.FixCorr <- function(evmod,Sims,SimID) {
  if (evmod$bifac_model) {
    fac_names <- names(evmod$model_struct$facs)
    phi.mat <- combinations(fac_names,k=2,replace=FALSE)
    phis <- apply(phi.mat,1,paste,collapse="")
  } else {
    phis <- evmod$model_struct$phis
  }
  
  rv <- list()
  for (phi in phis) {
    if (evmod$bifac_model) {
      value <- 0
    } else { # assumed 2 or more facs if this fn was called
      target <- sprintf("phiu_%s",phi)
      position <- match(SimID,Sims$sids)
      value <- pull(Sims,target)[position]
    }
    tmp.line <- sprintf("%s WITH %s@%f;"
                        ,str_sub(phi,1,1)
                        ,str_sub(phi,2,2)
                        ,value)
    rv <- c(rv,tmp.line)
  }
  rv <- c(rv," ")
  return(rv)
}

L.ModelEM.FixMeanVars <- function(evmod) {
  rv <- list()
  for (factor in names(evmod$model_struct$facs)) {
    lines.mv <- sprintf("[%s@0]; %s@1;",factor,factor)
    rv <- c(rv,lines.mv)
  }
  rv <- c(rv," ")
}

L.ModelEM.FixThresholds <- function(evmod,Sims,SimID) {
  rv <- list()
  for (dx in evmod$model_struct$useVars) {
    position <- match(SimID,Sims$sids)
    dxlab <- evmod$model_struct$varTypes$mylab[
      match(dx,evmod$model_struct$varTypes$oglab)
      ]
    colname <- sprintf("re_%s",dxlab)
    value <- pull(Sims,colname)[position]
    tmp.line <- sprintf("[%s$1@%f];",dx,value)
    rv <- c(rv,tmp.line)
  }
  rv <- c(rv," ")
  return(rv)
}

L.ModelEM.FixRegr <- function(evmod,SimID,EM_FullResults)  {
  mymod <- evmod$clone()
  y <- mymod$y
  xs <- mymod$xs
  p_id <- mymod$p_id
  p_uid <- sprintf("%s-%s",p_id,SimID)
  # Ok find now what xs should be fixed to.
  xs_fixed <- c()
  for (i in 1:length(xs)) {
    xvar <- xs[i]
    tgt_colname <- sprintf("ou_%s",xvar)
    tgt_row_index <- match(p_uid,EM_FullResults$uid)
    fix_val <- pull(EM_FullResults,tgt_colname)[tgt_row_index]
    xs_fixed <- c(xs_fixed
                  ,sprintf("%s@%f",xvar,fix_val))
  }
  line <- VarsToSplitLines(c(y,"ON",xs_fixed)
                           ,add.semicolon = TRUE
                           ,add.star = FALSE)
  rv <- list(line,"")
  return(rv)
}

L.ModelEM.RegrFSE <- function(evmod,SimID,EM_FullResults) {
  rv <- list("MODEL:","")
  if (!is_tibble(EM_FullResults)) {
    regrline <- L.ModelEM.RegrFull(evmod)
  } else {
    regrline <- L.ModelEM.FixRegr(evmod,SimID,EM_FullResults) 
  }
  rv <- c(rv,regrline)
  return(rv)
}

L.ModelEM.FixSEM <- function(evmod,Sims,SimID,EM_FullResults) {
  
  # Initialize the command list.
  rv <- list("MODEL:","")
  
  # Fix loadings
  lines.load <- L.ModelEM.FixLoad(evmod,Sims,SimID)
  rv <- c(rv,lines.load)
  
  # Fix thresholds
  lines.thr <- L.ModelEM.FixThresholds(evmod,Sims,SimID)
  rv <- c(rv,lines.thr)
  
  # Fix factor correlations
  # (If we have more than one factor present, we should deal
  # with "WITH" statements)
  if (length(names(evmod$model_struct$facs)) > 1) {
    lines.with <- L.ModelEM.FixCorr(evmod,Sims,SimID)
    rv <- c(rv,lines.with)
  }
  
  # Fix factor means and variances!
  lines.facmeanvar <- L.ModelEM.FixMeanVars(evmod)
  rv <- c(rv,lines.facmeanvar)
  
  # Finally, specify the regression we want to run!
  if (!is_tibble(EM_FullResults)) {
    regrline <- L.ModelEM.RegrFull(evmod)
  } else {
    regrline <- L.ModelEM.FixRegr(evmod,SimID,EM_FullResults)
  }
  rv <- c(rv,regrline," ")
  
  return(rv)
}

L.ModelEM.CatVars <- function(evmod) {
  # We will filter out to ensure that all the vars
  # included in this statement are dependent vars
  # (i.e., observed vars for estimating a factor or the y-var of
  # interest)
  mymod <- evmod$clone()
  dv <- mymod$catvars()
  dv <- dv[!(dv %in% mymod$xs)]
  rv <- L.CatVars(dv)
  return(rv)
}

# This function will create statements from USEVARIABLES through ANALYSIS,
# which will be invariant for a given model.
AssemblePrototype.EM <- function(evmod,KEY_VARS,estimator) {
  
  mymod <- evmod$clone()
  
  rv <- list()
  
  # Add USEVARIABLES
  tmp <- L.UseVars(mymod$usevars())
  rv <- c(rv,tmp)
  
  # Add CATEGORICAL statement
  # Note: we only include dependent variable (ys) AND variables that
  # will be used in estimating the latent factors in this statement.
  tmp <- L.ModelEM.CatVars(mymod)
  rv <- c(rv,tmp)
  
  # Add ID/STRAT/CLUSTER/WT section
  sample <- mymod$sample
  tmp <- L.IDSpecs(sample,KEY_VARS)
  rv <- c(rv,tmp)
  
  # Add ANALYSIS section
  tmp <- L.Analysis(estimator)
  rv <- c(rv,tmp)
  
  return(rv)
  
}

AssembleFile.EM <- function(MplusPrefix,mod,KEY_VARS,Sims,SimID,estimator="WLSMV"
                            ,EM_FullResults=NA) {
  
  # Copy our model structure to avoid weirdness
  evmod <- mod$clone()
  
  # Create structure to hold lines
  rv <- list()
  
  # Mplus prefix contains NAMES and DATA statement.
  # We need to change the file path to reference prior folder.
  lines.prefix <- L.EditPrefixPath.EM(MplusPrefix)
  rv <- c(rv,lines.prefix)
  
  # Get USEVARIABLES through ANALYSIS by calling this function:
  tmp <- AssemblePrototype.EM(evmod,KEY_VARS,estimator)
  rv <- c(rv,tmp)
  
  # Add the MODEL section
  lat_method <- evmod$lat_method
  if ((is.na(lat_method)) | (lat_method == "FSE")) {
    tmp <- L.ModelEM.RegrFSE(evmod,SimID,EM_FullResults)
    rv <- c(rv,tmp)
  } else if (lat_method == "FIXED_SEM") {
    tmp <- L.ModelEM.FixSEM(evmod,Sims,SimID,EM_FullResults)
    rv <- c(rv,tmp)
  }
  
  # Add OUTPUT section
  tmp <- "OUTPUT: SAMPSTAT CINT STANDARDIZED;"
  rv <- c(rv,tmp,"")
  
  rv <- unlist(rv)
  names(rv) <- NULL
  
  return(rv)
}


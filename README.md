# zero-cases-pfactor

Files in this repository correspond to our manuscript in prep, *Spurious empirical support for the p-factor arises with the inclusion of undiagnosed cases.* 

Preprint available here: https://psyarxiv.com/4tazx/

OSF page: https://osf.io/c9zys/

## Requirements

- Mplus version 8 or above installed on your machine and accessible using the MplusAutomation package from your R environment
- R version 3.6.3

### Important notice about R environment and reproducibility

This code is meant to be run using R 3.6.3. An `renv` environment is included in the repository for reproducibility, indicating all package dependencies. You will need to ensure that you are running the scripts from within the appropriate R environment and using the appropriate version, or alternatively, manually reconfigure your default environment to match the description specified in our `renv` container. For more information and tutorials on `renv`, consult the following page: https://rstudio.github.io/renv/articles/renv.html

## Code description

The code included here is used to generate a folder with all the necessary files to explore how dropping undiagnosed subjects impacts parameter estimates, for one sample, for one model. The file you will need to leverage is `ScriptGenerator.R`. Edit this file to ensure you correctly specify:

- `ScriptFolder`: point to the location of this repository
- `TargetFolderRoot`: where you would like the folder with scripts you'll need to be located
- `RunningFolderRoot`: where the analyses will be ran, if different (e.g., path in your computer cluster)
- `SampleID`: what sample you're looking at
- `AnalysisType`: whether you will just do CFA, CFA+EFA, or EFA
- `ModelID`: the model specification you would like to examine (see what we included)
- Other options, including dropping undiagnosed cases based on a specific diagnosis, etc.

This means that, if you wanted to look at how dropping zero cases impacted NESARC Wave 1 loadings in a correlated three factor model.

Once the folder is generated, you would run the R Script in the folder. If you are using RStudio and have loaded correctly the `renv` environment, you can simplify things by relying on the `renv::run()` function and calling the script directly from RStudio as a job.

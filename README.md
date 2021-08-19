# zero-cases-pfactor

Files in this repository correspond to our manuscript in prep, *Spurious empirical support for the p-factor arises with the inclusion of undiagnosed cases.* 

Preprint available here: https://psyarxiv.com/4tazx/

OSF page: https://osf.io/c9zys/

## Code description

This code is meant to be run using R 3.6. An `renv` environment is included in the repository for reproducibility, indicating all package dependencies. 

The code included here is used to generate a folder with all the necessary files to explore how dropping undiagnosed subjects impacts parameter estimates, for one sample, for one model. The file you will need to leverage is `ScriptGenerator.R`. Edit this file to ensure you correctly specify:

- `ScriptFolder`: point to the location of this repository
- `TargetFolderRoot`: where you would like the folder with scripts you'll need to be located
- `RunningFolderRoot`: where the analyses will be ran, if different (e.g., path in your computer cluster)
- `SampleID`: what sample you're looking at
- `AnalysisType`: whether you will just do CFA, CFA+EFA, or EFA
- `ModelID`: the model specification you would like to examine (see what we included)
- Other options, including dropping undiagnosed cases based on a specific diagnosis, etc.

This means that, if you wanted to look at how dropping zero cases impacted NESARC Wave 1 loadings in a correlated three factor model.

Once the folder is generated, you would run the R Script in the folder.
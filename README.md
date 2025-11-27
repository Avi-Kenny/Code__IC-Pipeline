# Controlled Risk Curves 

## Reproducibility

General setup:

- Download the repository from https://github.com/Avi-Kenny/Code__IC-Pipeline/

- Assume that we have a high performance computing environment with a slurm scheduler.

- Assume that we have R 4.4.2 loaded, e.g. on the Fred Hutch cluster, ml R-bundle-CRAN/2024.11-foss-2024a, which includes gsl

- Load additional modules
  ```{bash}
  ml CMake/3.29.3-GCCcore-13.3.0
  ml GLPK/5.0-GCCcore-13.3.0
  ml Pandoc/2.13
  ml ICU/75.1-GCCcore-13.3.0
  ```

- Assume that we have renv installed. If not, open R console at the project level (the folder containing this readme file), and run the following commands to install the current renv from CRAN, which is 1.1.5 as of September 2025. A different renv version may also work.
  ```{r}
  install.packages("renv")
  
  packageVersion("renv")  
  ```

- Run the following R command at the project level to install package dependencies. Note that superlearner and pch need to be installed separately as shown below.
  ```{R}
  renv::restore() # choose restore when presented with options
   ```



### ENSEMBLE trial severe correlates manuscript

Make sure that the analysis string in R/analysis.R is set to "Janssen (partA)" as shown below.
```{r}
  cfg2 <- list(analysis="Janssen (partA)", calc_ests=T, seed=1)
```

Modify R/analysis.R so that cfg2$folder_local below _if (cfg2$analysis=="Janssen (partA)") {_ points to the local copy of analysis-ready data file, e.g.,
  ```{r}
    cfg2$folder_local <- "../covpn/adata/"
  ```

To generate the controlled risk plots, run the following commands in a bash shell on the repo root level. 64 jobs are created to run in parallel on a slurm cluster. Each job renders one of the 64 set of plots and supporting files. The script run_r.sh contains the following code:
```{bash}
sbatch --array=1-64 run_r.sh
```

When all jobs are done, the controlled risk and VE plots will be in the folder Figures + Tables\Janssen (partA) plots. To get the mediation results, run the following command in R at the project level:
```{r}
source("make_mediation_table.R")
```



### Sanofi stage 2 correlates manuscript

Make sure that the analysis string in R/analysis.R is set to "Sanofi" as shown below.
```{r}
  cfg2 <- list(analysis="Sanofi", calc_ests=T, seed=1)
```

Modify R/analysis.R so that cfg2$folder_local below _if (cfg2$analysis=="Sanofi") {_ points to the local copy of analysis-ready data file, e.g.,
  ```{r}
    cfg2$folder_local <- "../covpn/adata/"
  ```

To generate the controlled risk plots, run the following commands in a bash shell on the repo root level. 90 jobs are created to run in parallel on a slurm cluster. Each job renders one of the 90 set of plots and supporting files. The script run_r.sh contains the following code:
```{bash}
sbatch --array=1-90 run_r.sh
```


# Controlled Risk Curves 

## Reproducibility

General setup:

- Download code from https://github.com/Avi-Kenny/Code__IC-Pipeline/. Different projects may have different releases.

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

Use this release: https://github.com/Avi-Kenny/Code__IC-Pipeline/archive/refs/tags/v1.0.zip

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
Mapping of the output:
```
  TID	 Arm	  Mrk	TS	Type	Mth
 1-15	Vacc	 1-15	 1	Abs		6
16-30	Vacc	16-30	 1	OverB	6
31-45	Vacc	31-45	 2	Abs		5
46-60	Vacc	46-60	 2	OverB	5
61-75	Plac	 1-15	 1	Abs		6
76-90	Plac	31-45	 2	Abs		5
91-105	Vacc	31-45	 2	Abs		5
106-120	Vacc	46-60	 2	OverB	5

  TID	 Arm	  Mrk	TS	Type	M
 1-15	Vacc	 1-15	 1	Abs		6	V+P set 1
61-75	Plac	 1-15	 1	Abs		6	V+P set 1
31-45	Vacc	31-45	 2	Abs		5	V+P set 2 ***** Changed from M6 to M5 (2024-10-22)
76-90	Plac	31-45	 2	Abs		5	V+P set 2 ***** Changed from M6 to M5 (2024-10-22)



PDF MERGES

TID
 1-30	2024-08-13 Controlled Risk TS1 M6 (Cox + NP)
 1-30	2024-08-13 Controlled Vaccine Efficacy TS1 M6 (Cox + NP)
31-60	2024-08-13 Controlled Risk TS2 M6 (Cox + NP)
31-60	2024-08-13 Controlled Vaccine Efficacy TS2 M6 (Cox + NP)
91-120	2024-08-13 Controlled Risk TS2 M5 (Cox + NP)
91-120	2024-08-13 Controlled Vaccine Efficacy TS2 M5 (Cox + NP)
1*-15*	2024-08-13 Controlled Risk V+P TS1 M6 (NP)
16*-30*	2025-04-18 Controlled Risk V+P TS2 M5 (NP)



PAPER FIGURE NUMBERS

Plot		Fig		Pair	id_v	id_p		Y-max
V+P 8		5 A		1		8		68			6.0%
V+P 12		5 B		1		12		72			6.0%
i 53		5 C		2		53					3.5%
i 57		5 D		2		57					3.5%
V+P 1		S28 A	3		1		61			6.0%
V+P 2		S28 B	3		2		62			6.0%
i 46 		S28 C	4		46					6.5%
i 47		S28 D	4		47					6.5%
V+P 10		S29 A	5		10		70			6.0%
V+P 11		S29 B	5		11		71			6.0%
i 55		S29 C	6		55					3.0%
```
i 56		S29 D	6		56					3.0%

1111111

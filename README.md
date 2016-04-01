# panel-dpm2016

This code was downloaded from < https://github.com/compops/panel-dpm2016 > or from < http://liu.johandahlin.com/ > and contains the code used to produce the results in the technical report

* J. Dahlin, R. Kohn and T. B. Schön, *Bayesian inference for mixed effects models with heterogeneity*. Technical report LiTH-ISY-R-3091, Department of Electrical Engineering, Linköping University. March 2016. 

The papers are available as a preprint from < http://liu.johandahlin.com/ >.

Requirements
--------------
The code is written and tested for R 3.2.3. To run the code, you need to have mvtnorm, lme4 and RColorBrewer installed, which amounts to executing the command "install.packages(c("mvtnorm","RColorBrewer","lme4"))". 

Main script files
--------------
Please, make sure that the working directory is the same as the directory in which the following .R-files are found. Otherwise, the code will not run as the subroutines implementing the Gibbs samplers cannot be loaded.

**example1-mixture-normals/example1-run-comparison.R** This script recreates the example in Section 4.1 together with Figure 1.

**example2-mixedeffects/example2-run-comparison.R** This script estimates the model in Section 4.2 for a specific N and T. It was used for verfying the code before running the simulation study. Hence, the settings can differ slightly from the ones indicated in the report. 

**example2-mixedeffects/example2-run-gamma.R** This script recreates the simulation study presented in Section 4.2. It is run from the command line by executing "R CMD BATCH --no-save --no-restore '--args jj=XXX' example2-run-gamma.R" for XXX from 1 to 40 to get 40 independent runs. Please, make sure that the working directory set in the file is correct before executing the script. Also, make sure to run the data generating script (discussed below) before running this file.

**example3-reactiontest/example3-reactiontest.R** This script replicats the results of Setion 4.3.

**example2-mixedeffects/example2-generate-datasets.R** This script generates the data sets for the simulation study in Section 4.2.

**example2-mixedeffects/example2-computeLL.R** This script computes Tables 1 from the output of the simulation script (discussed above).

**example2-mixedeffects/example2-mixedeffects-plots.R** This script recreates Figure 2 from the output of the simulation script (discussed above).

**example3-reactiontest/example3-data.R** This script recreates Figure 3.

Supporting files
--------------
**example1-mixture-normals/samplers-normalmixture.R** Gibbs samplers for the finite mixture model and the DPM model in Section 4.1.

**example2-mixedeffects/samplers-mixedeffects-1dim.R** Gibbs samplers for the finite mixture model and the DPM model in Section 4.2.

**example3-reactiontest/samplers-mixedeffects-sleepmodel.R** Gibbs samplers for the finite mixture model and the DPM model in Section 4.3.


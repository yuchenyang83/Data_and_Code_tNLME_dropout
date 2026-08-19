# tNLME models with non-ignorable dropout
Supplement: Data and Code for "Modeling HIV Viral Dynamics Using a Nonlinear Mixed-Effects Framework for Heavy-Tailed Data with Informative Dropout" by Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang*

# Author responsible for the code #
For questions, comments or remarks about the code please contact responsible author, Yu-Chen Yang (d107053002@smail.nchu.edu.tw).

# Configurations #
The code was written/evaluated in R with the following software versions:
- R version 4.6.0 (64-bit)
- Platform: x86_64-w64-mingw32/x64 (64-bit)
- Running under: Windows 11 x64 (build 26200)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.utf8  LC_CTYPE=Chinese (Traditional)_Taiwan.utf8    LC_MONETARY=Chinese (Traditional)_Taiwan.utf8
[4] LC_NUMERIC=C                                  LC_TIME=Chinese (Traditional)_Taiwan.utf8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] nlme_3.1-169  mvtnorm_1.4-2 tidyr_1.3.2   scales_1.4.0  ggrepel_0.9.8 dplyr_1.2.1   ggplot2_4.0.3

loaded via a namespace (and not attached):
 [1] vctrs_0.7.3        cli_3.6.6          rlang_1.3.0        purrr_1.2.2        generics_0.1.4     S7_0.2.2           glue_1.8.1         labeling_0.4.3
 [9] tibble_3.3.1       lifecycle_1.0.5    compiler_4.6.0     RColorBrewer_1.1-3 Rcpp_1.1.2         pkgconfig_2.0.3    farver_2.1.2       lattice_0.22-9
[17] R6_2.6.1           tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.5     tools_4.6.0        withr_3.0.3        gtable_0.3.6

# Descriptions of the codes #
Please extract **Data_and_Code.zip** so that `master.R` and the subfolders `./Code`, `./function`, `./Data`, and `./Result` are located in the same current working directory.
The `getwd()` function returns the absolute pathname of the current working directory, which in our setup is typically set to `D:/Data_and_Code/`.

Before running the codes, install the following R packages:

    install.packages("mvtnorm")  # Version: 1.4-2
    install.packages("ggplot2")  # Version: 4.0.3
    install.packages("dplyr")    # Version: 1.2.1
    install.packages("ggrepel")  # Version: 0.9.8
    install.packages("scales")   # Version: 1.4.0
    install.packages("tidyr")    # Version: 1.3.2
    install.packages("nlme")     # Version: 3.1-169

R codes for the implementation of our methodology for an ACTG 398 dataset and simulation results are provided.

## File: ./master.R ##
**master.R** is the main reproduction script. It loads the required packages and sequentially sources the scripts that reproduce Figures 1, 2, B.1, B.2, B.3, and E.1 and Tables 1, 2, B.1, B.2, B.3, D.1, D.2, and F.1 from the stored fitted objects and intermediate simulation results.

The master script does not rerun the computationally intensive ACTG 398 model fitting or Monte Carlo simulation studies.

# Subfolder: ./function #
`./function` contains the functions used for fitting the NLME and tNLME models and for supporting the ACTG 398 analyses:

- (1) **tNLMMmissingSAEM.R**: main SAEM function for fitting the tNLME model with incomplete longitudinal responses under MCAR, MAR, or MNAR.
- (2) **NLMMmissingSAEM.R**: main SAEM function for fitting the NLME model with incomplete longitudinal responses under MCAR, MAR, or MNAR.
- (3) **analyze_realdata_AIDS.R**: supporting functions for the ACTG 398 analysis, including the nonlinear mean model, derivatives, correlation matrices, likelihood calculations, missing-response updates, and related computational steps.
- (4) **tNLMMmissingSAEM_sensitivity.R**: tNLME SAEM function used for the Metropolis-Hastings tuning sensitivity analysis reported in Appendix F.
- (5) **NLMMmissingSAEM_sensitivity.R**: NLME SAEM function used for the corresponding Metropolis-Hastings tuning sensitivity analysis.
- (6) **analyze_realdata_AIDS_sensitivity.R**: supporting functions used by the Appendix F tuning-sensitivity analysis.
- (7) **multiplot.R**: utility function for arranging multiple `ggplot2` objects in a common layout.

## Subfolder: ./function/fix_alpha ##
`./function/fix_alpha` contains the functions used for the fixed-alpha2 sensitivity analysis:

- (1) **tNLMMmissingSAEM_fixalpha.R**: modified tNLME SAEM function for fitting the MNAR model with alpha2 fixed at a prespecified value.
- (2) **run_sensitivity_analysis.R**: supporting functions for the fixed-alpha2 sensitivity analysis reported in Figure 2.

## Subfolder: ./function/simulation ##
`./function/simulation` contains the functions used in the simulation studies:

- (1) **tNLMMmissingSAEM.R**: main SAEM function for fitting the tNLME model in the simulation studies.
- (2) **NLMMmissingSAEM.R**: main SAEM function for fitting the NLME model in the simulation studies.
- (3) **simulate_dropout_study.R**: supporting functions for generating nonlinear longitudinal data with dropout and for constructing quantities required by the simulation-fitting routines.

# Subfolder: ./Code #
`./Code` contains all model-fitting, simulation, diagnostic, figure, and table scripts:

- (1) **fit_actg398.R**: main computational script for fitting the NLME and tNLME models under MCAR, MAR, and MNAR across random-effects Scenarios (I)-(III) and four within-subject correlation structures: UNC, CS, AR(1), and MA(1). The script also contains the Scenario (III) fits based on pooled nonlinear least-squares (`nls`) initial values for Appendix D and the 5000-iteration Scenario (III) AR(1) tNLME MNAR fit used in Appendix E.
- (2) **fit_actg398_fixalpha.R**: main computational script for the Figure 2 sensitivity analysis. It fits the selected Scenario (III) AR(1) tNLME model under MNAR with alpha2 fixed over the prespecified grid from -6 to 6.
- (3) **fit_actg398_sensitivity.R**: main computational script for Appendix F. It refits the Scenario (III) AR(1) tNLME model under MCAR, MAR, and MNAR for the Metropolis-Hastings tuning settings reported in Table F.1.
- (4) **simSAEM25.R**: main script for the primary heavy-tailed simulation with approximately 25% dropout across sample sizes N = 25, 50, 100, 200, and 400.
- (5) **simSAEM50.R**: main script for the primary heavy-tailed simulation with approximately 50% dropout across sample sizes N = 25, 50, 100, 200, and 400.
- (6) **simSAEM75.R**: main script for the primary heavy-tailed simulation with approximately 75% dropout across sample sizes N = 25, 50, 100, 200, and 400.
- (7) **simSAEMt25-50.R**: main script for the additional nearly Gaussian simulation (nu = 100) with approximately 25% dropout across sample sizes N = 25, 50, 100, 200, and 400.
- (8) **simSAEMt50-50.R**: main script for the additional nearly Gaussian simulation (nu = 100) with approximately 50% dropout across sample sizes N = 25, 50, 100, 200, and 400.
- (9) **simSAEMt75-50.R**: main script for the additional nearly Gaussian simulation (nu = 100) with approximately 75% dropout across sample sizes N = 25, 50, 100, 200, and 400.
- (10) **diagnostic_information.R**: general diagnostic utility for an already fitted NLME or tNLME object. It returns the approximated observed-data log-likelihood trajectory, parameter trace plots, parameter-estimate summaries, Metropolis-Hastings acceptance-rate summaries, and the mapping between the reported parameters and the stored parameter-trace columns.
- (11) **Fig1.R**: main script for reproducing Figure 1.
- (12) **Fig2.R**: main script for reproducing Figure 2.
- (13) **FigB1.R**: main script for reproducing Figure B.1.
- (14) **FigB2.R**: main script for reproducing Figure B.2.
- (15) **FigB3.R**: main script for reproducing Figure B.3.
- (16) **FigE1.R**: main script for reproducing Figure E.1, which displays the 5000-iteration approximated observed-data log-likelihood trajectory used to assess the stopping criterion.
- (17) **Tab1.R**: main script for reproducing Table 1.
- (18) **Tab2.R**: main script for reproducing Table 2.
- (19) **TabB1.R**: main script for reproducing Table B.1.
- (20) **TabB2.R**: main script for reproducing Table B.2.
- (21) **TabB3.R**: main script for reproducing Table B.3.
- (22) **TabD1.R**: main script for reproducing Table D.1 under the alternative pooled-`nls` initialization.
- (23) **TabD2.R**: main script for reproducing Table D.2 under the alternative pooled-`nls` initialization.
- (24) **TabF1.R**: main script for reproducing Table F.1, which summarizes the Metropolis-Hastings tuning sensitivity analysis.

### Diagnostic information ###
The diagnostic function can be applied to any compatible fitted object after the corresponding `.RData` file has been loaded. For example:

```r
load("D:/Data_and_Code/Data/fit_III_result.RData")
source("D:/Data_and_Code/Code/diagnostic_information.R")

gg = diagnostic_information(fit.t.III.ARp.MNAR)

gg$loglik.plot
gg$trace.plot
gg$diagnostic.table
gg$parameter.table
gg$acceptance.table
gg$acceptance.by.component
gg$trace.column.map
```

To inspect a Scenario (I) or Scenario (II) fit, first load **fit_I_result.RData** or **fit_II_result.RData**, respectively.

# Subfolder: ./Data #
`./Data` contains the source data, stored fitted objects, sensitivity-analysis results, and simulation results required by the reproduction scripts.

The files stored directly in `./Data` are:

- (1) **fit_I_result.RData**: collects the fitted NLME and tNLME models for random-effects Scenario (I) across the missingness mechanisms and within-subject correlation structures used in Table 1.
- (2) **fit_II_result.RData**: collects the fitted NLME and tNLME models for random-effects Scenario (II) across the missingness mechanisms and within-subject correlation structures used in Table 1.
- (3) **fit_III_result.RData**: collects the fitted models for random-effects Scenario (III), including the objects used in Figure 1 and Tables 1-2.
- (4) **fit_III_result_nls.RData**: collects the Scenario (III) NLME and tNLME fits obtained using the alternative pooled-`nls` initialization for Tables D.1 and D.2.
- (5) **fit.t.III.ARp.MNAR5000.RData**: stores the selected Scenario (III) AR(1) tNLME MNAR fit run for up to 5000 SAEM iterations for Appendix E and Figure E.1.

The subfolders in `./Data` are described below.

## Subfolder: ./Data/source ##
`./Data/source` contains:

- (1) **actg398.txt**: the ACTG 398 longitudinal HIV-1 RNA dataset used in the real-data analysis. The file contains 2,626 longitudinal observations from 481 participants and the variables `patid`, `nnrti`, `calwk`, `txday`, `logrna`, `cens`, and `trtarm`.

## Subfolder: ./Data/fix_alpha ##
`./Data/fix_alpha` contains the fitted Scenario (III) AR(1) tNLME MNAR objects used for the fixed-alpha2 sensitivity analysis in Figure 2.

The fitted-object filenames have the form

`fit.t.III.ARp.MNAR_alpha[value].RData`,

where `[value]` is one of the following fixed alpha2 values:

`-6, -4, -2, -1, -0.5, -0.1, -0.05, -0.01, -0.001, -0.0001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 4, 6`.

The folder also contains:

- **fixed_alpha.txt**: intermediate summary containing the fixed alpha2 values, fixed-effect estimates, standard errors, and confidence-limit information used by **Fig2.R**.

## Subfolder: ./Data/sensitivity_SAEM ##
`./Data/sensitivity_SAEM` contains the stored Scenario (III) AR(1) tNLME fits used for the Metropolis-Hastings tuning sensitivity analysis in Appendix F.

The fitted-object filenames have the form

`fit.t.III.ARp.[mechanism]_k[K]_c[c].RData`,

where:

- `[mechanism]` is `MCAR`, `MAR`, or `MNAR`;
- `[K]` is `5`, `10`, or `15`;
- `[c]` is `1.8`, `2.4`, or `3`.

Thus, the current sensitivity analysis contains 27 fitted objects: 3 missingness mechanisms x 3 values of K x 3 values of c.

The folder also contains:

- **TableF1_raw.txt**: intermediate file containing the extracted K value, missingness mechanism, proposal scaling constant c, number of estimated parameters, approximated observed-data log-likelihood, and overall Metropolis-Hastings acceptance rate used to construct Table F.1.

The implementation details for these tuning settings are provided in **fit_actg398_sensitivity.R**.

## Subfolder: ./Data/simulation ##
`./Data/simulation` contains six simulation-result folders:

- (1) **SS-simulation-t25**: primary heavy-tailed simulation results with approximately 25% dropout.
- (2) **SS-simulation-t50**: primary heavy-tailed simulation results with approximately 50% dropout.
- (3) **SS-simulation-t75**: primary heavy-tailed simulation results with approximately 75% dropout.
- (4) **SS-simulation-t25-50**: additional nearly Gaussian simulation results (nu = 100) with approximately 25% dropout.
- (5) **SS-simulation-t50-50**: additional nearly Gaussian simulation results (nu = 100) with approximately 50% dropout.
- (6) **SS-simulation-t75-50**: additional nearly Gaussian simulation results (nu = 100) with approximately 75% dropout.

Each of these six folders contains five subfolders, **SIM1**, **SIM2**, **SIM3**, **SIM4**, and **SIM5**, corresponding to sample sizes N = 25, 50, 100, 200, and 400, respectively.

### Primary simulation folders ###
The subfolders **SS-simulation-t25**, **SS-simulation-t50**, and **SS-simulation-t75** contain the simulation results for the primary heavy-tailed setting.

Each `SIM1`-`SIM5` subfolder contains:

- **aic.txt**: per-replication AIC values.
- **bic.txt**: per-replication BIC values.
- **drop_rate.txt**: realized dropout proportions.
- **iter.txt**: stored iteration-count information.
- **loglik.txt**: per-replication approximated maximized observed-data log-likelihood values.
- **MSE.txt**: per-replication mean squared prediction-error summaries.
- **MAE.txt**: per-replication mean absolute prediction-error summaries.
- **MAPE.txt**: per-replication mean absolute percentage prediction-error summaries.
- **MAR.para.est.txt**, **MCAR.para.est.txt**, and **MNAR.para.est.txt**: tNLME parameter estimates under the corresponding missingness specification.
- **MAR.para.fit.txt**, **MCAR.para.fit.txt**, and **MNAR.para.fit.txt**: NLME parameter estimates under the corresponding missingness specification.
- **MAR.se.est.txt**, **MCAR.se.est.txt**, and **MNAR.se.est.txt**: information-matrix standard errors from the corresponding tNLME fits.
- **MAR.se.fit.txt**, **MCAR.se.fit.txt**, and **MNAR.se.fit.txt**: information-matrix standard errors from the corresponding NLME fits.
- **realpara.txt**: true parameter values used to generate the simulated datasets.

### Additional nearly Gaussian simulation folders ###
The subfolders **SS-simulation-t25-50**, **SS-simulation-t50-50**, and **SS-simulation-t75-50** contain the additional nearly Gaussian simulation results with nu = 100.

Each `SIM1`-`SIM5` subfolder contains:

- **aic.txt**: per-replication AIC values.
- **bic.txt**: per-replication BIC values.
- **drop_rate.txt**: realized dropout proportions.
- **iter.txt**: stored iteration-count information.
- **loglik.txt**: per-replication approximated maximized observed-data log-likelihood values.
- **MSE.txt**: per-replication mean squared prediction-error summaries.
- **MAE.txt**: per-replication mean absolute prediction-error summaries.
- **MAPE.txt**: per-replication mean absolute percentage prediction-error summaries.
- **MNAR.iter.lnL.est.txt**: stored iteration-wise approximated observed-data log-likelihood trajectories from the tNLME fits under MNAR.
- **MNAR.iter.lnL.fit.txt**: stored iteration-wise approximated observed-data log-likelihood trajectories from the NLME fits under MNAR.
- **MNAR.para.est.txt**: tNLME parameter estimates under MNAR.
- **MNAR.para.fit.txt**: NLME parameter estimates under MNAR.
- **MNAR.se.est.txt**: information-matrix standard errors from the tNLME fits under MNAR.
- **MNAR.se.fit.txt**: information-matrix standard errors from the NLME fits under MNAR.
- **realpara.txt**: true parameter values used to generate the simulated datasets.
- **seed.txt**: stored random-seed and replication information.

# Subfolder: ./Result #
`./Result` contains all figures and tables generated by the reproduction scripts:

- (1) **Figure1.pdf**: displays visit-specific distributions of log10(RNA), treatment-group mean trajectories, visit-level outliers, subject-level outlier trajectories identified from the selected model, and the numbers and proportions of observed and missing responses at each scheduled visit.
- (2) **Figure2.pdf**: displays the fixed-effect estimates and 95% confidence intervals from the selected Scenario (III) AR(1) tNLME MNAR model across the fixed-alpha2 sensitivity grid.
- (3) **FigureB1.pdf**: shows one representative simulated dataset with N = 100, including subject-specific trajectories, response boxplots, and observed-response counts and percentages by treatment group.
- (4) **FigureB2.pdf**: shows the MSEs of the estimated parameters from the fitted tNLME model under MNAR across three dropout rates and sample sizes N = 25, 50, 100, 200, and 400.
- (5) **FigureB3.pdf**: compares the MSEs of parameter estimates from the NLME and tNLME models under MNAR in the additional nearly Gaussian simulation.
- (6) **FigureE1.pdf**: displays the approximated observed-data log-likelihood trajectory for the selected Scenario (III) AR(1) tNLME MNAR model over 5000 SAEM iterations.
- (7) **Table1.csv**: summarizes the 72 candidate NLME and tNLME models using the number of model parameters, maximized approximated observed-data log-likelihood, AIC, and BIC across Scenarios (I)-(III), four within-subject correlation structures, and three missingness mechanisms.
- (8) **Table2.csv**: reports ML estimates, standard errors, and absolute estimate-to-SE ratios for the Scenario (III) AR(1) tNLME models under MCAR, MAR, and MNAR.
- (9) **TableB1.csv**: compares information-matrix standard errors (IM SE) with empirical Monte Carlo standard deviations (MC SD) from the fitted tNLME model under MNAR across sample sizes and dropout rates.
- (10) **TableB2.csv**: reports average AIC, BIC, and MSPE values together with selection frequencies for the NLME and tNLME models under MNAR across sample sizes and dropout rates.
- (11) **TableB3.csv**: reports average AIC and BIC values together with selection frequencies for the tNLME model fitted under MCAR, MAR, and MNAR across sample sizes and dropout rates.
- (12) **TableD1.csv**: summarizes model-selection results for the Scenario (III) models fitted using the alternative pooled-`nls` initialization.
- (13) **TableD2.csv**: reports parameter estimates, standard errors, and absolute estimate-to-SE ratios for the Scenario (III) AR(1) tNLME fits obtained using the alternative pooled-`nls` initialization.
- (14) **TableF1.csv**: reports the approximated observed-data log-likelihood, AIC, BIC, and overall Metropolis-Hastings acceptance rate across the reported K and c tuning settings.

# Additional Remarks #
- Note (1): **master.R** reproduces the reported figures and tables from the stored fitted objects and intermediate simulation files. The computationally intensive fitting and simulation scripts are not sourced by `master.R`.
- Note (2): **fit_actg398.R**, **fit_actg398_fixalpha.R**, and **fit_actg398_sensitivity.R** can require substantial computation because they refit the ACTG 398 models. The fitted objects needed for reproducing the reported figures and tables are already stored under `./Data`.
- Note (3): **simSAEM25.R**, **simSAEM50.R**, **simSAEM75.R**, **simSAEMt25-50.R**, **simSAEMt50-50.R**, and **simSAEMt75-50.R** can also require substantial computation. The intermediate numerical results used by Figures B.1-B.3 and Tables B.1-B.3 are already stored under `./Data/simulation`.
- Note (4): To diagnose an individual fitted model, load the `.RData` file containing that model, source **diagnostic_information.R**, call `diagnostic_information(fit_object)`, and then display the desired returned plot or table component.
- Note (5): All paths in the reproduction workflow are defined relative to the repository root returned by `getwd()`. In our setup, the repository root is `D:/Data_and_Code/`.

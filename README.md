# tNLME models with non-ignorable dropout
Supplement: Data and Code for "Robust HIV Viral Dynamics: A Nonlinear Mixed-Effects Framework for Heavy-Tailed Data with Informative Dropout" by Yu-Chen Yang, Tsung-I Lin, Luis M. Castro, and Wan-Lun Wang*

# Author responsible for the code #
For questions, comments or remarks about the code please contact responsible author, Yu-Chen Yang (d107053002@smail.nchu.edu.tw).

# Configurations #
The code was written/evaluated in R with the following software versions:
   - R version 4.3.1 (2023-06-16 ucrt)
   - Platform: x86_64-w64-mingw32/x64 (64-bit)
   - Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.utf8  LC_CTYPE=Chinese (Traditional)_Taiwan.utf8    LC_MONETARY=Chinese (Traditional)_Taiwan.utf8
[4] LC_NUMERIC=C                                  LC_TIME=Chinese (Traditional)_Taiwan.utf8  

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] glue_1.7.0    pROC_1.18.5   gridExtra_2.3 ggtext_0.1.2  rlang_1.1.4   cowplot_1.1.3 mvtnorm_1.3-1 ggplot2_3.5.1 nlme_3.1-164

loaded via a namespace (and not attached):
 [1] crayon_1.5.3      vctrs_0.6.5       cli_3.6.3         generics_0.1.3    labeling_0.4.3    colorspace_2.1-1  plyr_1.8.9       
 [8] gridtext_0.1.5    scales_1.3.0      fansi_1.0.6       munsell_0.5.1     tibble_3.2.1      lifecycle_1.0.4   compiler_4.4.1   
[15] dplyr_1.1.4       Rcpp_1.0.13       pkgconfig_2.0.3   rstudioapi_0.16.0 farver_2.1.2      lattice_0.22-6    R6_2.5.1         
[22] tidyselect_1.2.1  utf8_1.2.4        pillar_1.9.0      magrittr_2.0.3    tools_4.4.1       withr_3.0.1       gtable_0.3.5     
[29] xml2_1.3.6 

# Descriptions of the codes # 
Please extract the file "Data_and_Code_tNLME_dropout.zip" to the "current working directory" of the R package.
The getwd() function returns the absolute pathname of the "current working directory", which in our setup is typically set to D:/Data_and_Code/.

Before running all of the codes, one needs to install the following R packages:
    
    install.packages("mvtnorm")  # Version: 1.3-1
    install.packages("ggplot2") # Version: 3.5.1
    install.packages("gridExtra") # Version: 2.3
    install.packages("ggtext") # Version: 0.1.2
    install.packages("rlang") # Version: 1.1.4
    install.packages("cowplot") # Version: 1.1.3
    install.packages("glue") # Version: 1.7.0
    install.packages("nlme") # Version: 3.1-164
    

R codes for the implementation of our methodology for an ACTG 398 dataset and simulation results are provided.

# Data and Code for AIDS data #
## Subfolder: ./Data/source ##
`./Data/source`
contains
- **actg398.txt**: the dataset from the ACTG 398 study.


## Subfolder: ./function ##
`./function`
contains the program (function) of
       
- (1) **tNLMMmissingSAEM.R** main script for fitting the tNLME model with non-ignorable dropout;
- (2) **NLMMmissingSAEM.R** main script for fitting the NLME model with non-ignorable dropout;
- (3) **multiplot.R** main script for combining multiple plots by ggplot2 package;
- (4) **analyze_realdata_AIDSR.R** that can assist in fitting the NLME and tNLME models with non-ignorable dropout.

## Subfolder: ./function/fixed_alpha ##
`./function/fixed_alpha`
contains the program (function) of
	
- (1) **tNLMMmissingSAEM_fixalpha.R** main script for fitting the tNLME model with MNAR dropout;
- (2) **run_sensitivity_analysis.R** that assists in fitting the tNLME model with non-ignorable dropout

## Subfolder: ./Code ##
`./Code`
contains the program (function) of
- (1) **fit_actg398.R** is the main script for fitting the NLME and tNLME models under the three mechanisms. It includes four structures for within-patient autocorrelation.
- (2) **fit_actg398_fixalpha.R** is the main script for fitting the tNLME model under the MNAR mechanism, where the parameter alpha2 is fixed at values ranging from -6 to 6.
- (3) **Fig1.R** main script for reproducing Figure 1;
- (4) **Fig2.R** main script for reproducing Figure 2; 
- (5) **Tab1.R** main script for Table 1;
- (6) **Tab2.R** main script for Table 2;

## Subfolder: ./Data ##
`./Data`
contains
- (1) **fit_result.RData** collects the fitting results of the NLME and tNLME models with the three missing data mechanisms and four within-patient autocorrelation structures applied to the ACTG 398 dataset.
- (2) **fixed_alpha.txt** contains the results of the sensitivity analysis, collecting fixed alpha 2 values within the range of -6, -4, -2, -1, -0.5, -0.1, -0.05, -0.01, -0.001, -0.0001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 4, 6.

## Subfolder: ./Result ##
`./Result`
contains
- (1) **Figure1.eps** Figure 1 displays, in the top panel, the longitudinal trajectories of $\log_10(RNA)$ responses from 481 patients, overlaid with boxplots at each scheduled visit (Day 0, 2, 4, 8, 16, and 24). The bottom panel presents the corresponding bar chart of the number and percentage of observed responses at each time point, illustrating the extent of missing data across the study period.
- (2) **Figure2.eps** displays estimated fixed effects parameters (solid line) and 95% confidence intervals (dashed lines) from the fitted tNLME model with MNAR.
- (3) **Table1.csv** lists the fitting results for the 24 candidate models, including the number of unknown parameters, maximized log-likelihood values together with AIC and BIC scores for determining the preferred model.
- (4) **Table2.csv** presents a comparison of the ML estimates of parameters together with their standard errors (SE) obtained from the tNLME model with AR(1) errors under the three mechanisms.


# Data and Code for Simulation #
## Subfolder: ./function/simulation ##
`./function/simulation`
contains the program (function) of
- (1) **tNLMMmissingSAEM.R** main script for fitting the tNLME model with non-ignorable dropout;
- (2) **NLMMmissingSAEM.R** main script for fitting the NLME model with non-ignorable dropout;
- (3) **simulate_dropout_study.R** that can assist for fitting the NLME and tNLME models with non-ignorable dropout

## Subfolder: ./Code ##
`./Code`
contains the program (function) of
- (1) The main script **simSAEM25.R** is used to fit NLME and tNLME models under three missing data mechanisms to simulated datasets with a 25% dropout rate.
- (2) The main script **simSAEM50.R** is used to fit NLME and tNLME models under three missing data mechanisms to simulated datasets with a 50% dropout rate.
- (3) The main script **simSAEM75.R** is used to fit NLME and tNLME models under three missing data mechanisms to simulated datasets with a 75% dropout rate.
- (4) **FigB1.R** main script for reproducing Figure B.1;
- (5) **FigB2.R** main script for reproducing Figure B.2; 
- (6) **TabB1.R** main script for Table B.1;
- (7) **TabB2.R** main script for Table B.2;
- (8) **TabB3.R** main script for Table B.3;


## Subfolder: ./Data/simulation ##
`./Data/simulation`
	subfolder contains
- (1) The subfolder `SS-simulation-t25` contains the fitting results of the NLME and tNLME models for data with a 25% dropout rate, where three missing data mechanisms are based on 100 repetitions across various sample sizes.
- (2) The subfolder `SS-simulation-t50` contains the fitting results of the NLME and tNLME models for data with a 50% dropout rate, where three missing data mechanisms are based on 100 repetitions across various sample sizes.
- (3) The subfolder `SS-simulation-t75` contains the fitting results of the NLME and tNLME models for data with a 75% dropout rate, where three missing data mechanisms are based on 100 repetitions across various sample sizes.

## Subfolder: ./Data/simulation/SS-simulation-t25 ##
`./Data/simulation/SS-simulation-t25`
	subfolder contains

- (1) 5 subsubfolders: `./SIM1`, `./SIM2`, `./SIM3`, and `./SIM4`, and `./SIM5`,  storing the **aic.txt**, **bic.txt**, **drop_rate.txt**, **iter.txt**, **loglik.txt**, **MSE.txt**, **MAE.txt**, **MAPE.txt**, **MAR.para.est.txt**, **MAR.para.fit.txt**, **MAR.se.est.txt**, **MAR.se.fit.txt**, **MCAR.para.est.txt**, **MCAR.para.fit.txt**, **MCAR.se.est.txt**, **MCAR.se.fit.txt**, **MNAR.para.est.txt**, **MNAR.para.fit.txt**, **MNAR.se.est.txt**, **MNAR.se.fit.txt**, and **realpara.txt** text files for Appendix B. Simulations;

## Subfolder: ./Data/simulation/SS-simulation-t50 ##
`./Data/simulation/SS-simulation-t50`
	subfolder contains
	
- (1) 5 subsubfolders: `./SIM1`, `./SIM2`, `./SIM3`, and `./SIM4`, and `./SIM5`,  storing the **aic.txt**, **bic.txt**, **drop_rate.txt**, **iter.txt**, **loglik.txt**, **MSE.txt**, **MAE.txt**, **MAPE.txt**, **MAR.para.est.txt**, **MAR.para.fit.txt**, **MAR.se.est.txt**, **MAR.se.fit.txt**, **MCAR.para.est.txt**, **MCAR.para.fit.txt**, **MCAR.se.est.txt**, **MCAR.se.fit.txt**, **MNAR.para.est.txt**, **MNAR.para.fit.txt**, **MNAR.se.est.txt**, **MNAR.se.fit.txt**, and **realpara.txt** text files for Appendix B. Simulations;

## Subfolder: ./Data/simulation/SS-simulation-t75 ##
 `./Data/simulation/SS-simulation-t75`
	subfolder contains
	
- (1) 5 subsubfolders: `./SIM1`, `./SIM2`, `./SIM3`, and `./SIM4`, and `./SIM5`,  storing the **aic.txt**, **bic.txt**, **drop_rate.txt**, **iter.txt**, **loglik.txt**, **MSE.txt**, **MAE.txt**, **MAPE.txt**, **MAR.para.est.txt**, **MAR.para.fit.txt**, **MAR.se.est.txt**, **MAR.se.fit.txt**, **MCAR.para.est.txt**, **MCAR.para.fit.txt**, **MCAR.se.est.txt**, **MCAR.se.fit.txt**, **MNAR.para.est.txt**, **MNAR.para.fit.txt**, **MNAR.se.est.txt**, **MNAR.se.fit.txt**, and **realpara.txt** text files for Appendix B. Simulations;

## Subfolder: ./Result ##
`./Result`
contains
- (1) **FigureB1.eps** shows the trajectory plot of responses for one simulated case of size N = 100 (left panel), box plots for responses for subjects (right top panel), and bar chart of frequency distribution of number of observed responses  in the two prespecified groups (right bottom panel).
- (2) **FigureB2.eps** shows the MSE scores for the estimated parameters under fitted tNLME model with the MNAR mechanism for three dropout rates across various sample sizes.
- (3) **TableB1.csv** lists the simulation results for assessing the asymptotic standard errors (IM SE) and empirical standard deviations (MC Sd) of parameters estimates under fitted tNLME model with the MNAR mechanism across various sample sizes.
- (4) **TableB2.csv** reports the average AIC, BIC, and MSPE scores together with frequencies (in parentheses) supported by the two criteria for the NLME and tNLME models with various sample sizes and dropout rates. The result for the best performance per row is highlighted in bold.
- (5) **TableB3.csv** summarizes the average AIC and BIC values, together with the frequencies (in parentheses) of model selection based on these criteria, for the tNLME models fitted under the MCAR, MAR, and MNAR mechanisms across different sample sizes and dropout proportions.


## Additional Remark ##
- Note (1): One can directly run each "source(.)" described in **master.R** file in the separate R session to obtain the results.
- Note (2): The fitting results of the considered models for the ACTG 398 dataset obtained by running **fit_actg398.R** and **fit_actg398_fixalpha.R** have been stored in **fit_result.RData** and **fixed_alpha.txt**.
- Note (3): To reproduce Tables 1-2 and Figure 1-2 in Section 4, please load the 'fit_result.RData' and 'fixed_alpha.txt' file in subfolder `./Data/`, and then run **Fig1-2.r** and **Tab1-2.r** script in subfolder `./Code/`.
- Note (4): Since **fit_actg398_fixalpha.R** and **fit_actg398.R** takes a long time to run, to reproduce Tables 1-2 and Figure 1-2 in Section 4, please load the **fit_result.RData** and **fixed_alpha.txt** file in subfolder `./Data/`, and then run **Fig1-2.r** and **Tab1-2.r** script in subfolder `./Code/`.
- Note (5): Since **simSAEM25.R**, **simSAEM50.R**, and **simSAEM75.R** takes a long time to run, to reproduce numerical results in Appendix B. Simulations, we record these intermediate numerical results so that one can use the R codes **FigB1.r**, **FigB2.r**, **TabB1.r**, **TabB2.r**, and **TabB3.r** to obtain the final results based on files stored in `./simulation/SS-simulation-t25/SIM1`, `./simulation/SS-simulation-t25/SIM2`, `./simulation/SS-simulation-t25/SIM3`, `./simulation/SS-simulation-t25/SIM4`, `./simulation/SS-simulation-t25/SIM5`, `./simulation/SS-simulation-t50/SIM1`, `./simulation/SS-simulation-t50/SIM2`, `./simulation/SS-simulation-t50/SIM3`, `./simulation/SS-simulation-t50/SIM4`, `./simulation/SS-simulation-t50/SIM5`, `./simulation/SS-simulation-t75/SIM1`, `./simulation/SS-simulation-t75/SIM2`, `./simulation/SS-simulation-t75/SIM3`, `./simulation/SS-simulation-t75/SIM4`, and `./simulation/SS-simulation-t75/SIM5` subfolders.

Anton Nakov, "Jackknife Instrumental Variables Estimation: 
Replication and Extension of Angrist, Imbens and Krueger (1999)",
Journal of Applied Econometrics

The zip file contains nine Ox programs:

1. moncarlo.ox 

Runs the Monte Carlo simulations for the five models. The parameters are set to replicate the original results of Angrist, Imbens and Krueger (1999).

To replicate the new results, set the sample size "T" to 25; the autocorrelation "rho" to 0.95; or change the covariance matrix "var" for Model 3.

2. schooling.ox

Estimates returns to schooling. The program is set to replicate the original results of Angrist, Imbens and Krueger (1999). Set "samples" to 100 to replicate the new results.

3. ols.ox

Performs ordinary least squares estimation. Input: dependent variable y, matrix of regressors X. Output: point estimate and standard error of the slope coefficients 


4. tsls.ox

Performs two-stage least squares estimation. Input: dependent variable y, matrix of endogenous regressors X, matrix of instruments Z. Output: point estimate and standard error of the slope coefficients in the main regression.


5. liml.ox

Performs limited information maximum likelihood estimation. Input: dependent variable y1, matrix of endogenous regressors Y2, matrix of included exogenous regressors X, matrix of excluded exogenous regressors Z. Output: point estimate and standard error of the slope coefficients in the main regression.

6. jive1.ox

Implements the JIVE1 estimator of Angrist, Imbens and Krueger (1999).
Input: dependent variable y, matrix of endogenous regressors X, matrix of instruments Z. Output: point estimate and standard error of the slope coefficients in the main regression.

7. jive2.ox

Implements the JIVE2 estimator of Angrist, Imbens and Krueger (1999).
Input: dependent variable y, matrix of endogenous regressors X, matrix of instruments Z. Output: point estimate and standard error of the slope coefficients in the main regression.

8. stats.ox

Calculates the quantiles, Median Absolute Error, and 95% confidence interval coverage


9. setseed.ox

A function by Charles Bos to reset the iSeed. Input: integer iSeed, value for ranseed if iSeed > 0, else current time is used. Output: none, ranseed is reset

In addition, there are six Excel files with output tables:

table I.xls  : replication of Table I (Monte Carlo simulations) in Angrist, Imbens and Krueger (1999).

table II.xls : replication of Table II (Returns to schooling) in Angrist, Imbens and Krueger (1999).

table 1.xls  : Table 1 of the replication paper containing Anton Nakov's simulations of Model 3 (where the discrepancy with the original paper was found)

table 2.xls  : Table 2 of the replication paper: Returns to schooling estimated from 100 random sub-samples of 50,000 observations 

table 3.xls  : Table 3 of the replication paper:  Monte Carlo simulations with a small sample size, N = 25

table 4.xls  : Table 4 of the replication paper: Autocorrelated errors, rho = 0.95


The data used is the same as that in the original article and can be downloaded from the JAE Data Archive: 

http://qed.econ.queensu.ca/jae/

Angrist, J, Imbens, G,  Krueger A. 1999. Jackknife instrumental variables estimation. Journal of Applied Econometrics 14: 57-67. 



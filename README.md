# BDDN: Bayesian Dynamic Differential Network

This repository contains the implementation of the **Bayesian Dynamic Differential Network (BDDN)** analysis pipeline, which models time-varying protein interactions under drug treatment using a covariance regression framework with Bayesian inference.

## 📁 Repository Structure

- `code/covreg.cpp`:  
  C++ implementation of the Gibbs sampling algorithm for the covariance regression model.

- `code/covreg.R`:  
  R wrapper to call `covreg.cpp` via Rcpp and perform posterior sampling for parameters **B** and **Ψ**.  
  **Output**: Posterior samples of `B` (regression coefficients) and `Ψ` (noise covariance).

- `code/prereg.R`:  
  Post-processing script to compute:
  - **Ω(d,t)**: precision matrix,
  - **R(d,t)**: partial correlation matrix,
  - **D(t)**: differential score matrix,

  based on the output of `covreg.R`.

- `code/BDDN.R`:  
  Final analysis pipeline including:
  - Functional clustering of edge trajectories,
  - Construction and visualization of dynamic differential networks.

## 📊 Example Data

Example datasets are provided in the `data/` folder:

- `data/Y.csv`: **Protein expression matrix**  
  - Rows = samples, Columns = proteins  
  - This is the main response matrix used in the model.

- `data/XdatS.csv`: **Sample metadata**  
  - Columns:  
    - `UID`: unique sample identifier  
    - `treatment`: 0 (control) or 1 (treated)  
    - `hours`: observation time point  
  - This is used as covariates in the regression model.

You can load the data as:

```r
XdatS = read.csv("data/XdatS.csv", stringsAsFactors = FALSE)
Y = as.matrix(read.csv("data/Y.csv", row.names = 1, check.names = FALSE))

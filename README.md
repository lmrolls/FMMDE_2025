# Macroeconomic Forecasting with FMMDE
Reproducibility Package for Macroeconomic Forecasting using Factor Models with Martingale Difference Errors

## Package Information
- **Date Assembled**: April 23, 2025
- **Author(s)**: Luca Mattia Rolla, Alessandro Giovannelli
- **Contact**: Luca Mattia Rolla (lmrolla92@gmail.com, University of Rome “Tor Vergata”, Italy), Alessandro Giovannelli (alessandro.giovannelli@univaq.it, University of L'Aquila, Italy)
- **Repository URL**: https://github.com/lmrolls/Macroeconomic-Forecasting-with-FMMDE

## Overview
This repository contains the code and documentation to reproduce the results in "Macroeconomic Forecasting using Factor Models with Martingale Difference Errors" by Luca Mattia Rolla and Alessandro Giovannelli, submitted to *International Journal of Forecasting*. The package reproduces the paper's tables (1–5) and figures.

## Repository Structure
- `code/`: MATLAB scripts and functions.
  - `codeForSimulations/`: Scripts for simulations.
    - `DGPs/`: Data-generating process functions.
      - `fLeeShao.m`: Generates simulated data for a three-dimensional linear factor model as described in Example 5 of Lee & Shao (2018).
      - `fLeeShaoNonLinear.m`: Generates simulated data for a three-dimensional nonlinear factor model as described in Example 6 of Lee & Shao (2018).
    - `factorEstimation/`: Factor estimation methods.
      - `factorLAM.m`: Estimates a factor model using the cumulative linear covariance matrix with spectral decomposition and eigenvalue ratio method, as described in Lam et al. (2011) and Lam & Yao (2012).
      - `factorMDDM.m`: Estimates a factor model using the cumulative Martingale Difference Divergence Matrix (CMDDM) with spectral decomposition and eigenvalue ratio method, as described in Lee & Shao (2018).
      - `factorMDDM4.m`: Variant of MDDM method.
      - `factorMDDMtab4.m`: MDDM method specific to Table 4.
      - `fixedFactorLAM.m`: Estimates a factor model with a fixed number of factors using the cumulative linear covariance matrix, as described in Lam et al. (2011) and Lam & Yao (2012).
      - `fixedFactorsMDDM.m`: Estimates a factor model with a fixed number of factors using the cumulative Martingale Difference Divergence Matrix (CMDDM), as described in Lee & Shao (2018).
      - `MDDM.m`: Computes the Martingale Difference Divergence Matrix (MDDM) between two vector time series, as described in Lee & Shao (2018).
      - `seqTest.m`: Implements a sequential testing procedure to determine the number of latent factors in FMMDE models, based on methods in Lee & Shao (2018), Wang et al. (2022), and Lam & Yao (2011).
      - `stockwatson2002.m`: Factor estimation based on Stock and Watson (2002).
    - `Table4/`: Scripts for Table 4.
      - `MainFile_TAB4.m`: Main script.
      - `Tab4.xlsx`: Output file.
      - `subfunctions/standardize.m`: Helper function for data standardization.
    - `Table5/`: Scripts for Table 5.
      - `MainFile_TAB5.m`: Main script.
      - `Table5_MSFE_Ratios.xlsx`: Output file.
      - `subfunctions/standardize.m`: Helper function for data standardization.
    - `Tables123/`: Scripts for Tables 1–3.
      - **Description**: Contains scripts to generate Tables 1–3 via Monte Carlo simulations evaluating factor model estimation accuracy under varying parameters (true number of factors, data-generating process, signal-to-noise ratio, and sample sizes). The main script conducts simulations with nested loops over sample sizes (n=50, T=50 or T=100), true number of factors (r=2 to 5), DGPs (1 to 3), signal-to-noise parameters (theta: 0.5*r, r, 3*r, 5*r), and 1000 Monte Carlo replications using `SimulFun.m`. It employs Sequential Testing and Eigenvalue Ratio methods, with parameters like MDDM (k0=1) and p-value (0.05). Results include empirical probabilities displayed in the console (matching LaTeX table format) and exported to Excel.
      - `MainFile_TAB1.asv`: Auto-saved version of Table 1 script.
      - `MainFile_TAB1.m`: Script specific to Table 1.
      - `MainFile_TAB123.m`: Main script for Tables 1–3 Monte Carlo simulations.
      - `format.xlsx`: Output file for table formatting.
      - `Simulation_Results.xlsx`: Output file containing simulation results for both sample sizes (n=50, T=50 and n=50, T=100).
      - `subfunctions/factorSim.m`: Helper function for factor simulations.
      - `subfunctions/SimulFun.m`: Core simulation function.
      - `subfunctions/standardize.m`: Helper function for data standardization.
  - `plotsLuca/`: Scripts and data for figures.
    - `nfactors.mat`: Data for number of factors plot.
    - `nfactorsBAI.mat`: Data for BAI factors plot.
    - `nfactorsEIG.mat`: Data for eigenvalue-based factors plot.
    - `plots.m`: Script to generate figures.
- `README.md`: This file.
- `LICENSE`: MIT License.

## Computing Environment for Simulations
- **Software**: MATLAB R2023b or later
- **Toolboxes**:
  - Statistics and Machine Learning Toolbox (12.4)
  - Linear Algebra (MATLAB core)
  - Parallel Computing Toolbox
- **License**: MIT License (see `LICENSE`)
- **Hardware**: Tested on a desktop computer with the following specifications:
  - **CPU**: Intel Core i7-9700 (8 cores, 8 threads, 3.00 GHz base, up to 4.70 GHz turbo)
  - **RAM**: 16 GB DDR4-2666
  - **Storage**: 512 GB SSD
  - **Operating System**: Windows 10
- **Setup Instructions**:
  1. Install MATLAB R2021a or later.
  2. Ensure required toolboxes are installed (see above).
  3. Clone this repository: `git clone https://github.com/lmrolls/Macroeconomic-Forecasting-with-FMMDE.git` or download as a ZIP.
  4. Open MATLAB and set the working directory to `code/`:
     ```matlab
     cd('path/to/Macroeconomic-Forecasting-with-FMMDE/code');
     ```
  5. Add paths manually or create a `startup.m` script (not included) to add all subfolders:
     ```matlab
     addpath(genpath('codeForSimulations'));
     ```
## Computing Environment for the Empirical Analysis

## Data
- **Sharable Data**: Not explicitly included in the provided structure. If datasets (e.g., `.mat` files) are required, they should be placed in a `data/` folder with descriptions.
- **Non-Sharable Data**: If applicable, include access instructions in a `data/external/non_sharable_data.md` file.
- **Intermediary Data**: Outputs like `Tab4.xlsx`, `Table5_MSFE_Ratios.xlsx`, and `Simulation_Results.xlsx` are generated by the scripts and saved in their respective folders.

## Running the Reproducibility Check
- **Computer Used**: Desktop computer with Intel Core i7-9700 (8 cores, 8 threads, 3.00 GHz base, up to 4.70 GHz turbo), 16 GB DDR4-2666 RAM, 512 GB SSD, Windows 10.
- **Instructions**:
  1. Clone the repository: `git clone https://github.com/lmrolls/Macroeconomic-Forecasting-with-FMMDE.git`.
  2. Open MATLAB and set the working directory to `code/`:
     ```matlab
     cd('path/to/Macroeconomic-Forecasting-with-FMMDE/code');
     ```
  3. Add paths:
     ```matlab
     addpath(genpath('codeForSimulations'));
     ```
  4. Run analysis scripts:
     - For Tables 1–3:
       ```matlab
       run('codeForSimulations/Tables123/MainFile_TAB123.m');
       ```
       Outputs: `codeForSimulations/Tables123/Simulation_Results.xlsx`
     - For Table 4:
       ```matlab
       run('codeForSimulations/Table4/MainFile_TAB4.m');
       ```
       Output: `codeForSimulations/Table4/Tab4.xlsx`
     - For Table 5:
       ```matlab
       run('codeForSimulations/Table5/MainFile_TAB5.m');
       ```
       Output: `codeForSimulations/Table5/Table5_MSFE_Ratios.xlsx`
     - For Figures:
       ```matlab
       run('plotsLuca/plots.m');
       ```
       Outputs: Figures generated from `nfactors.mat`, `nfactorsBAI.mat`, `nfactorsEIG.mat`
  5. Outputs are saved in their respective folders (e.g., `Tables123/`, `Table4/`, `Table5/`).

## Mapping Code to Paper Outputs
- **Tables 1–3**: Generated by `codeForSimulations/Tables123/MainFile_TAB123.m` using functions like `factorSim.m` and `SimulFun.m`. Output: `Simulation_Results.xlsx`.
- **Table 4**: Generated by `codeForSimulations/Table4/MainFile_TAB4.m` using `standardize.m`. Output: `Tab4.xlsx`.
- **Table 5**: Generated by `codeForSimulations/Table5/MainFile_TAB5.m` using `standardize.m`. Output: `Table5_MSFE_Ratios.xlsx`.
- **Figures**: Generated by `plotsLuca/plots.m` using `nfactors.mat`, `nfactorsBAI.mat`, and `nfactorsEIG.mat`.

## Notes
- Ensure MATLAB and required toolboxes are installed.
- Contact lmrolla92@gmail.com or alessandro.giovannelli@univaq.it for questions or issues.

## References
- Lam, C., & Yao, Q. (2012). Factor modeling for high-dimensional time series: Inference for the number of factors. *The Annals of Statistics*, 40(2), 694–726. DOI: 10.1214/12-AOS970
- Lam, C., Yao, Q., & Bathia, N. (2011). Estimation of latent factors for high-dimensional time series. *Biometrika*, 98(4), 901–918. DOI: 10.1093/biomet/asr047
- Lee, C. E., & Shao, X. (2018). Martingale Difference Divergence Matrix and its application to dimension reduction for stationary multivariate time series. *Journal of the American Statistical Association*, 113(521), 216–229. DOI: 10.1080/01621459.2016.1240082
- Wang, G., Zhu, K., & Shao, X. (2022). Testing for the Martingale Difference Hypothesis in Multivariate Time Series Models.

## License
This package is licensed under the MIT License (see `LICENSE`).

## Acknowledgments
The authors are grateful to Tommaso Proietti for his valuable and insightful comments, which led to several improvements in both the presentation and the content of the paper. The authors are also grateful to the participants of the 41st International Symposium on Forecasting 2022. Alessandro Giovannelli gratefully acknowledges financial support from the Italian Ministry of Education, University and Research, Progetti di Ricerca di Interesse Nazionale, research project 2020-2023, project 2020N9YFFE.
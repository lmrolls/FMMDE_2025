# Macroeconomic Forecasting with FMMDE
Reproducibility Package for Macroeconomic Forecasting using Factor Models with Martingale Difference Errors

## Package Information
- **Date Assembled**: April 23, 2025
- **Author(s)**: Luca Mattia Rolla, Alessandro Giovannelli
- **Contact**: Luca Mattia Rolla (lmrolla92@gmail.com, University of Rome “Tor Vergata”, Italy), Alessandro Giovannelli (alessandro.giovannelli@univaq.it, University of L'Aquila, Italy)
- **Repository URL**: https://github.com/lmrolls/Macroeconomic-Forecasting-with-FMMDE

## Overview
This repository contains the code and documentation to reproduce the results in "Macroeconomic Forecasting using Factor Models with Martingale Difference Errors" by Luca Mattia Rolla and Alessandro Giovannelli, submitted to *International Journal of Forecasting*. The package reproduces the paper's tables (1–5), which are based on simulations, and figures, which are part of the empirical analysis. The empirical analysis is forthcoming and not yet fully included.

## Repository Structure
- `code/`: MATLAB scripts and functions.
  - `codeForSimulations/`: Scripts for simulations (Tables 1–5).
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
      - `stockwatson2002.m`: Estimates a factor model using principal components, as described in Stock & Watson (2002a, 2002b).
    - `subfunctionsTAB123/`: Helper functions for Tables 1–3.
      - `factorSim.m`: Helper function for factor simulations.
      - `SimulFun.m`: Core simulation function.
      - `standardize.m`: Helper function for data standardization.
    - `Table1/`: Scripts for Table 1.
      - `MainFile_TAB1.m`: Script for Table 1 Monte Carlo simulations.
    - `Table2/`: Scripts for Table 2.
      - `MainFile_TAB2.m`: Script for Table 2 Monte Carlo simulations.
    - `Table3/`: Scripts for Table 3.
      - `MainFile_TAB3.m`: Script for Table 3 Monte Carlo simulations.
    - `Table4/`: Scripts for Table 4.
      - `MainFile_TAB4.m`: Main script.
      - `Tab4.xlsx`: Output file.
      - `subfunctions/standardize.m`: Helper function for data standardization.
    - `Table5/`: Scripts for Table 5.
      - `MainFile_TAB5.m`: Main script.
      - `Table5_MSFE_Ratios.xlsx`: Output file.
      - `subfunctions/standardize.m`: Helper function for data standardization.
  - `empiricalAnalysis/`: Scripts and data for empirical analysis (figures).
    - `plotsLuca/`: Scripts and data for figures.
      - `nfactors.mat`: Data for number of factors plot.
      - `nfactorsBAI.mat`: Data for BAI factors plot.
      - `nfactorsEIG.mat`: Data for eigenvalue-based factors plot.
      - `plots.m`: Script to generate figures.
- `outputFromSimulations/`: Output files from simulations.
  - `tables/`: Output tables.
    - `Tab_1.xlsx`: Output file containing simulation results for Table 1.
- `README.md`: This file.
- `LICENSE`: MIT License.

## Computing Environment for Simulations
- **Software**: MATLAB R2023b
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
  1. Install MATLAB R2023b.
  2. Ensure required toolboxes are installed (see above).
  3. Clone this repository: `git clone https://github.com/lmrolls/Macroeconomic-Forecasting-with-FMMDE.git` or download as a ZIP.
  4. Open MATLAB and set the working directory to `code/`:
     ```matlab
     cd('path/to/Macroeconomic-Forecasting-with-FMMDE/code');
     ```
  5. Add paths manually or create a `startup.m` script (not included) to add all subfolders:
     ```matlab
     addpath(genpath('codeForSimulations'));
     addpath(genpath('empiricalAnalysis'));
     ```

## Computing Environment for the Empirical Analysis
- **Note**: The empirical analysis is not yet fully included in this package. The `empiricalAnalysis/plotsLuca/` folder contains scripts and data for generating figures, but a complete setup (e.g., data sources, additional toolboxes) is not provided. Details will be provided in a future update. Contact the authors for more information.

## Data
- **Sharable Data**: Tables 1–5 are based on simulations, and data are generated by scripts in `codeForSimulations/DGPs/` (e.g., `fLeeShao.m`, `fLeeShaoNonLinear.m`). The figures in `empiricalAnalysis/plotsLuca/` may require external datasets (e.g., FRED-MD, as described in McCracken & Ng, 2016), which are not included in the repository. Details on accessing these datasets will be provided in a future update.
- **Non-Sharable Data**: Not applicable for the current simulation-based tables. If empirical data are added in the future, access instructions will be provided in a `data/external/non_sharable_data.md` file.
- **Intermediary Data**: Outputs like `Tab_1.xlsx`, `Tab4.xlsx`, and `Table5_MSFE_Ratios.xlsx` are generated by the simulation scripts and saved in `outputFromSimulations/tables/` (for Table 1) and `codeForSimulations/Table4/` or `Table5/` (for Tables 4–5).

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
     addpath(genpath('empiricalAnalysis'));
     ```
  4. Run analysis scripts:
     - For Table 1:
       ```matlab
       run('codeForSimulations/Table1/MainFile_TAB1.m');
       ```
       Output: `outputFromSimulations/tables/Tab_1.xlsx`
     - For Table 2:
       ```matlab
       run('codeForSimulations/Table2/MainFile_TAB2.m');
       ```
       Output: To be saved in `outputFromSimulations/tables/` (e.g., `Tab_2.xlsx`, not explicitly listed in provided structure).
     - For Table 3:
       ```matlab
       run('codeForSimulations/Table3/MainFile_TAB3.m');
       ```
       Output: To be saved in `outputFromSimulations/tables/` (e.g., `Tab_3.xlsx`, not explicitly listed in provided structure).
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
       run('empiricalAnalysis/plotsLuca/plots.m');
       ```
       Outputs: Figures generated from `nfactors.mat`, `nfactorsBAI.mat`, `nfactorsEIG.mat`
  5. Outputs are saved in their respective folders (e.g., `outputFromSimulations/tables/` for Table 1, `codeForSimulations/Table4/` for Table 4, `codeForSimulations/Table5/` for Table 5).

## Mapping Code to Paper Outputs
- **Table 1**: Generated by `codeForSimulations/Table1/MainFile_TAB1.m` using functions in `subfunctionsTAB123/` (`factorSim.m`, `SimulFun.m`, `standardize.m`). Output: `outputFromSimulations/tables/Tab_1.xlsx`.
- **Table 2**: Generated by `codeForSimulations/Table2/MainFile_TAB2.m` using functions in `subfunctionsTAB123/`. Output: Expected in `outputFromSimulations/tables/` (e.g., `Tab_2.xlsx`).
- **Table 3**: Generated by `codeForSimulations/Table3/MainFile_TAB3.m` using functions in `subfunctionsTAB123/`. Output: Expected in `outputFromSimulations/tables/` (e.g., `Tab_3.xlsx`).
- **Table 4**: Generated by `codeForSimulations/Table4/MainFile_TAB4.m` using `subfunctions/standardize.m`. Output: `codeForSimulations/Table4/Tab4.xlsx`.
- **Table 5**: Generated by `codeForSimulations/Table5/MainFile_TAB5.m` using `subfunctions/standardize.m`. Output: `codeForSimulations/Table5/Table5_MSFE_Ratios.xlsx`.
- **Figures**: Generated by `empiricalAnalysis/plotsLuca/plots.m` using `nfactors.mat`, `nfactorsBAI.mat`, and `nfactorsEIG.mat`.

## Notes
- Ensure MATLAB and required toolboxes are installed.
- Tables 1–5 are based on simulations using data generated by scripts in `codeForSimulations/DGPs/`. The figures in `empiricalAnalysis/plotsLuca/` are part of the empirical analysis, which is not yet fully included and may require external datasets (e.g., FRED-MD, McCracken & Ng, 2016) in a future update.
- Outputs for Tables 2 and 3 are not explicitly listed in the provided structure but are expected to be saved in `outputFromSimulations/tables/`. Verify output file names (e.g., `Tab_2.xlsx`, `Tab_3.xlsx`) after running the scripts.
- Contact lmrolla92@gmail.com or alessandro.giovannelli@univaq.it for questions or issues.

## References
- Lam, C., & Yao, Q. (2012). Factor modeling for high-dimensional time series: Inference for the number of factors. *The Annals of Statistics*, 40(2), 694–726. DOI: 10.1214/12-AOS970
- Lam, C., Yao, Q., & Bathia, N. (2011). Estimation of latent factors for high-dimensional time series. *Biometrika*, 98(4), 901–918. DOI: 10.1093/biomet/asr047
- Lee, C. E., & Shao, X. (2018). Martingale Difference Divergence Matrix and its application to dimension reduction for stationary multivariate time series. *Journal of the American Statistical Association*, 113(521), 216–229. DOI: 10.1080/01621459.2016.1240082
- McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research. *Journal of Business & Economic Statistics*, 34(4), 574–589.
- Stock, J. H., & Watson, M. W. (2002a). Forecasting using principal components from a large number of predictors. *Journal of the American Statistical Association*, 97(460), 1167–1179.
- Stock, J. H., & Watson, M. W. (2002b). Macroeconomic forecasting using diffusion indexes. *Journal of Business & Economic Statistics*, 20(2), 147–162.
- Wang, G., Zhu, K., & Shao, X. (2022). Testing for the Martingale Difference Hypothesis in Multivariate Time Series Models.

## License
This package is licensed under the MIT License (see `LICENSE`).

## Acknowledgments
The authors are grateful to Tommaso Proietti for his valuable and insightful comments, which led to several improvements in both the presentation and the content of the paper. The authors are also grateful to the participants of the 41st International Symposium on Forecasting 2022. Alessandro Giovannelli gratefully acknowledges financial support from the Italian Ministry of Education, University and Research, Progetti di Ricerca di Interesse Nazionale, research project 2020-2023, project 2020N9YFFE.
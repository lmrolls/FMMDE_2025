%%#########################################################################################
% CODE USED TO GENERATE TABLE 5: Linear Factors
%
% Purpose: Evaluate the out-of-sample one-step-ahead predictive performance of three factor
% estimation methods: FMMDE (Factor Model Maximum Distance Estimator), SW (Stock-Watson 2002),
% and LYB (Lee-Yu-Bai) using simulated data from a linear factor model.
%
% Parameters:
% - NN_start: Cross-sectional dimensions (N = [100, 300, 500])
% - TT_start: Time-series dimensions (T = [200, 500, 1000])
% - rr: True number of factors (rr = 3)
% - nreps: Number of Monte Carlo replications (400)
%
% Methodology:
% - Simulated data X_t (N-dimensional) are generated from a 3-dimensional factor model
%   (Eq. 5 in the paper) using fLeeShao(T,N,rr).
% - Factors F^{m}_t (3-dimensional) are estimated from X_t using methods m = {FMMDE, SW, LYB}.
% - For each variable X_{t,k} (k=1,...,N), forecast X_{t+1,k} using:
%     X_{t+1,k} = beta*F^{m}_t + epsi_t, where beta is estimated via OLS.
% - Compute the mean squared forecast error (MSFE) for each method by averaging squared
%   forecast errors across N variables per replication.
% - Calculate MSFE ratios: MSFE^{FMMDE}/MSFE^{m} (m = SW, LYB).
% - Average the MSFE ratios over 400 Monte Carlo replications.
%
% Dependencies:
% - fLeeShao(): Generates data from the 3-factor model (Eq. 5).
% - stockwatson2002(): Estimates factors using Stock-Watson (2002).
% - factorMDDM2(): Estimates factors using FMMDE.
% - factorLAM(): Estimates factors using LYB.
% - standardize(): Standardizes the input data.
%
% Output: Matrices of averaged MSFE ratios (MSFE^{FMMDE}/MSFE^{LYB}, MSFE^{FMMDE}/MSFE^{SW})
% for each (N,T) combination.
%%#########################################################################################

close all; clear all; clc;

addpath('../factorEstimation');
addpath('../DGPs');
addpath('subfunctions');

% Set random number generator for replicability
rng(1, 'twister');

% Parameters
NN_start = [100 300 500]; % Cross-sectional dimensions
TT_start = [200 500 1000]; % Time-series dimensions
rr = 3; % Number of factors
nreps = 400; % Number of Monte Carlo replications
k0 = 1;

% Initialize cell arrays for MSE storage
mse_cla_sw = cell(length(NN_start), length(TT_start));
mse_cla_md = cell(length(NN_start), length(TT_start));
mse_cla_lam = cell(length(NN_start), length(TT_start));

for TT = 1:length(TT_start)
    T = TT_start(TT);
    disp(['T = ', num2str(T)]);
    for NN = 1:length(NN_start)
        N = NN_start(NN);
        disp(['N = ', num2str(N)]);
        
        % Pre-allocate arrays for forecasts and true values
        yfor_sw = zeros(nreps, N);
        yfor_md = zeros(nreps, N);
        yfor_lam = zeros(nreps, N);
        ytrue = zeros(nreps, N);
        
        parfor rep = 1:nreps
            disp(rep)
            substream = RandStream('mt19937ar', 'Seed', rep);
            RandStream.setGlobalStream(substream);

            % Disable rank-deficient matrix warning
            warning('off', 'MATLAB:rankDeficientMatrix');

            % FACTOR DATA GENERATION
            x = fLeeShao(T, N, rr);
            x = standardize(x);
            Xt = x(1:end-1, :); % X_t
            Xout = x(end, :); % X_{t+1}

            % FACTOR ESTIMATION
            [Fmddm] = factorMDDM2(Xt, k0, rr); % FMMDE
            [Fsw] = stockwatson2002(Xt, rr); % SW
            [Flam] = factorLAM(Xt, k0, rr); % LYB

            for k = 1:N
                yt = Xt(:, k);
                ytout = Xout(:, k);
                ysub = yt(2:end, 1);

                % Lagged Factors (F_{t-1})
                Fswsub = Fsw(1:end-1, :);
                Fmdsub = Fmddm(1:end-1, :);
                Flamsub = Flam(1:end-1, :);

                % Factors at time t (F_t)
                Fsfor = Fsw(end, :);
                Fmdfor = Fmddm(end, :);
                Flamfor = Flam(end, :);

                % OLS ESTIMATION
                b_hat_sw = Fswsub\ysub;
                yfor_sw(rep, k) = Fsfor*b_hat_sw;
                
                b_hat_md = Fmdsub\ysub;
                yfor_md(rep, k) = Fmdfor*b_hat_md;
                
                b_hat_lam = Flamsub\ysub;
                yfor_lam(rep, k) = Flamfor*b_hat_lam;
                
                ytrue(rep, k) = ytout; % Assign once
            end
        end

        % MSE COMPUTATION (outside parfor)
        mse_cla_sw{NN,TT} = mean((ytrue - yfor_sw).^2, 2); % Mean across N variables per replication
        mse_cla_md{NN,TT} = mean((ytrue - yfor_md).^2, 2); % Size: nreps x 1
        mse_cla_lam{NN,TT} = mean((ytrue - yfor_lam).^2, 2);
        
        % Average MSE across replications
        GM_mse_sw(NN,TT) = mean(mse_cla_sw{NN,TT});
        GM_mse_md(NN,TT) = mean(mse_cla_md{NN,TT});
        GM_mse_lam(NN,TT) = mean(mse_cla_lam{NN,TT});
    end
end

% Compute and display MSFE ratios
disp('MSFE^{FMMDE}/MSFE^{LYB}:');
GM_mse_md ./ GM_mse_lam
disp('MSFE^{FMMDE}/MSFE^{SW}:');
GM_mse_md ./ GM_mse_sw

% Save results
%save('table5_results.mat', 'GM_mse_sw', 'GM_mse_md', 'GM_mse_lam');
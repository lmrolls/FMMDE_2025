%%#################################################################################
% CODE USED TO GENERATE TABLE 5: Nonlinear Factors
%
% NN_start: different cross section dimensions to loop over
% TT_start: different time series dimensions to loop over
% k0_start: different values of k0 to loop over
%
% rR: TRUE NUMBER OF FACTORS,
% nreps: number of monte carlo replications
%
%
% MAKES USE OF THE FUNCTION fLeeShaoNonLinear() which generates data from
% the factor model in EQ.5
%  THE FUNCTION GENERATES DATA FROM THREE LATENT FACTORS
% 
%%#########################################################################################

close all; clear all; clc;

addpath('../factorEstimation');
addpath('../DGPs');
addpath('subfunctions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set random number generator for replicability
rng(1, 'twister');

% Number of Predictors
NN_start    = [100 300 500];
% Number of Observations
TT_start    = [200 500 1000];

k0=1;


for TT = 1:length(TT_start)
    T = TT_start(TT);
    disp(T)
    for NN = 1:length(NN_start)
        N = NN_start(NN);
        disp(N)
        parfor rep = 1:1000
            disp(rep)

            substream = RandStream('mt19937ar', 'Seed', rep);
            RandStream.setGlobalStream(substream);

            %=========================
            %  FACTOR DATA GENERATION
            %=========================

            x               = fLeeShaoNonLinear(T,N);    % generate NONLINEAR factors
            x               = standardize(x);
            Xt              = x(1:end-1,:);              % X_t
            Xout            = x(end,:);                  % X_{t+1} value to be forecasted

            %====================
            %  FACTOR ESTIMATION
            %====================

            [Fmddm]         = factorMDDM2(Xt, k0, rr);  % FMMDE factors
            [Fsw]           = stockwatson2002(Xt,rr);   % SW  factors
            [Flam]          = factorLAM2(Xt, k0, rr);   % LYB factors

            for k=1:N % Loop over the n simulated variables

                yt      = Xt(:,k);   % k-th lagged variable y_t = X_{t,k} becomes dependent variable yt
                ytout   = Xout(:,k); % y_{t+1} = X_{t+1,k} future variable to be forecasted

                ysub    = yt(2:end,1); % y_{t}

                % Lagged Factors for One Step Ahead Forecast: F_{t-1}
                Fswsub  = Fsw(1:end-1,:);
                Fmdsub  = Fmddm(1:end-1,:);
                Flamsub = Flam(1:end-1,:);

                % Factors at time t: F_t
                Fsfor   = Fsw(end,:);
                Fmdfor  = Fmddm(end,:);
                Flamfor = Flam(end,:);

                %==================
                %  OLS ESTIMATION
                %==================

                b_hat_sw         = Fswsub\ysub;    %        y_t = bhat * F_{t-1}
                yfor_sw(rep,k)   = Fsfor*b_hat_sw; % yhat_{t+1} = bhat * F_{t}
            
                b_hat_md         = Fmdsub\ysub;
                yfor_md(rep,k)   = Fmdfor*b_hat_md;
         
                b_hat_lam        = Flamsub\ysub;
                yfor_lam(rep,k)  = Flamfor*b_hat_lam;
                ytrue(rep,k)     = ytout;           % y_{t+1}
            end

        end

        %==================
        %  MSE COMPUTATION
        %==================


        GM_mse_sw(NN,TT)      = mean(mse_cla_sw{NN,TT});

        GM_mse_md(NN,TT)      = mean(mse_cla_md{NN,TT});

        GM_mse_lam(NN,TT)      = mean(mse_cla_lam{NN,TT});

    end
end


GM_mse_md./GM_mse_sw;
GM_mse_md./GM_mse_lam;
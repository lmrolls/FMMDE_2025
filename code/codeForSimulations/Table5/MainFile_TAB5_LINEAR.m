%%#################################################################################
% CODE USED TO GENERATE TABLE 5: Linear Factors
%
% NN_start: different cross section dimensions to loop over
% TT_start: different time series dimensions to loop over
%
% rr: TRUE NUMBER OF FACTORS,
% nreps: number of monte carlo replications
%
% This code evaluates the out of sample one step ahead predictive capabilities
% of the three methods SW, LYB and FMMDE
%
% The simulated data X_t, a N dimensional stochastic process,
% are obtained from the 3-dimensional factor process defined in equation 5 in th paper. 
% Be F^{m}_t the factors estimated from X_t using method m 
%
% We assume we know the true number of factor so that F^{m}_t is
% taken to be 3-dimensional.
%
% The forecasting equation takes this form y_{t+1} = beta*F^{m}_t +epsi_t
% where y_t = X_{t,n} so that we obtain N forecasts, one for each dimension
% in the stochastic process X_t. 
%
% For each method m we compute the mean one step ahead squared forecast
% error (MSFE)  averaging over the N forecast errors committed in the prediction of the N
% variables in X_{t}.
%
% The performance is evaluated observing the ratio between the MSFEs computed for different methods
% The MSFE of FMMDE is always at the numerator.
% We have therefore MSFE^{FMMDE}/MSFE^{m}, m being LYB or SW. 
%
% Finally, we average the ratio of the MSFEs over 1000 monte carlo
% replications to obtain the final output reported in the tables
%
% The code MAKES USE OF THE FUNCTION fLeeShao() which generates data from 
% the factor model in EQ.5 of the paper
%  THE FUNCTION GENERATES DATA FROM THREE LATENT FACTORS
% 
% Also uses the function stockwatson2002() 
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
rr = 3;

k0=1;

    for TT = 1:length(TT_start)
        T = TT_start(TT);
        disp(T)
        for NN = 1:length(NN_start)
            N = NN_start(NN);
            disp(N)
            parfor rep = 1:400
                disp(rep)

                substream = RandStream('mt19937ar', 'Seed', rep);
                RandStream.setGlobalStream(substream);

                %=========================
                %  FACTOR DATA GENERATION
                %=========================

                x               = fLeeShao(T,N,rr); % generate LINEAR factors
                x               = standardize(x);
                Xt              = x(1:end-1,:);     % X_t
                Xout            = x(end,:);         % X_{t+1} value to be forecasted

                %====================
                %  FACTOR ESTIMATION
                %====================

                [Fmddm]          = factorMDDM2(Xt, k0, rr); % FMMDE factors
                [Fsw]            = stockwatson2002(Xt,rr);  % SW  factors
                [Flam]           = factorLAM2(Xt, k0, rr);   % LYB factors

                for k=1:N % Loop over the n simulated variables

                    yt      = Xt(:,k);   % k-th lagged variable y_t = X_{t,k} becomes dependent variable yt
                    ytout   = Xout(:,k); % y_{t+1} = X_{t+1,k} future variable to be forecasted
      
                    ysub           = yt(2:end,1); % y_{t}

                    % Lagged Factors for One Step Ahead Forecast: F_{t-1}
                    Fswsub         = Fsw(1:end-1,:);
                    Fmdsub         = Fmddm(1:end-1,:);
                    Flamsub         = Flam(1:end-1,:);

                    % Factors at time t: F_t
                    Fsfor          = Fsw(end,:);
                    Fmdfor         = Fmddm(end,:);
                    Flamfor         = Flam(end,:);

                    %==================
                    %  OLS ESTIMATION
                    %==================
    
                    b_hat_sw        = Fswsub\ysub;  %        y_t = bhat * F_{t-1}
                    yfor_sw(rep,k)  = Fsfor*b_hat_sw;% yhat_{t+1} = bhat * F_{t}

                    b_hat_md        = Fmdsub\ysub;
                    yfor_md(rep,k)  = Fmdfor*b_hat_md;
                    ytrue(rep,k)    = ytout;
                    
                    b_hat_lam        = Flamsub\ysub;
                    yfor_lam(rep,k)  = Flamfor*b_hat_lam;
                    ytrue(rep,k)     = ytout;   % y_{t+1}
                end
              
            end

            %==================
            %  MSE COMPUTATION
            %==================
           
            mse_cla_sw{NN,TT}     = mean((ytrue - yfor_sw).^2);
            mse_cla_md{NN,TT}     = mean((ytrue  - yfor_md).^2);
            mse_cla_lam{NN,TT}     = mean((ytrue  - yfor_lam).^2);


            GM_mse_sw(NN,TT)      = mean(mse_cla_sw{NN,TT});

            GM_mse_md(NN,TT)      = mean(mse_cla_md{NN,TT});

            GM_mse_lam(NN,TT)     = mean(mse_cla_lam{NN,TT});

        end
    end


GM_mse_md./GM_mse_lam;
GM_mse_md./GM_mse_sw;
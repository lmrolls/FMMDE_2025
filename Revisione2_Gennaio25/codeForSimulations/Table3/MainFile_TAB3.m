%%#################################################################################
% CODE USED TO GENERATE TABLES FROM 1 TO 3 IN THE PAPER
% 
% THIS SCRIPT DETAILS THREE EXTERNAL LOOPS OVER THE SIMULATION PARAMETERS
% NAMELY: TRUE NUMBER OF FACTORS, DGP, THETA (SIGNAL TO NOISE PARAMETER)

% THE 4TH INTERNAL LOOP CONDUCTS THE MONTE CARLO SIMULATION'S REPLICATIONS
% THIS USES THE FUNCTION SuperSimulBari4 WHICH CONDUCTS DATA SIMULATION AND
% FACTOR ESTIMATION

% n:cross sectional dimension simulated data, T: time series dimension, k0: parameter MDDM
% pval: critical value used in sequential testing procedure
% mOutTest,mOutRAT: proportion of successful attempts in estimating of true number of factors
% Sequential Test and Eigenvalue Ratio methods respectively
%
% meanOutTest, meanOutRAT: mean of the estimated number of factors
% nIters: number of monte carlo exercise replications

% r: variable in the first loop, selects the true number of factors in the DGP
% DGP: Sets DGP from 1 to 3 in tables 1 to 3 in the paper 
% thetaMethod: the parameter theta in the paper taking values
% 1 ==> thetaMethod = 0.5*r
% 2 ==> thetaMethod = r
% 3 ==> thetaMethod = 3*r
% 4 ==> thetaMethod = 5*r

%%#########################################################################################
clear 
clc
addpath SW_MDDM
addpath numberF
addpath utils
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n      = 100;
T      = 100;
k0     = 1; 
pval   = 0.05; 


mOutTest = zeros(3,4,5) ;  
meanOutTest = zeros(3,4,5) ;

mOutRAT = zeros(3,4,5) ; 
meanOutRAT = zeros(3,4,5) ;

nIters = 1000;


for r = 2:5           % number of factors that are simulated
    for DGP = 1:3     % DGP from 1 to 3 
        
        for thetaMethod = 1:4   % Theta, the signal to noise parameter

            out = NaN(nIters,2);
            
            parfor k = 1:nIters  % number of Monte Carlo replications
                
                out(k,:) = SimulFun(T,n,k0,r,thetaMethod,DGP,pval);
                disp([k,thetaMethod,DGP,r]);
                
            end
            disp(mean(out==r))
           
            
            mOutTest(DGP,thetaMethod,r) = mean(out(:,1)==r);
            meanOutTest(DGP,thetaMethod,r) = mean(out(:,1));
            
            mOutRAT(DGP,thetaMethod,r) = mean(out(:,2)==r);
            meanOutRAT(DGP,thetaMethod,r) = mean(out(:,2)); 
        
            
        end
        
    end
    
end



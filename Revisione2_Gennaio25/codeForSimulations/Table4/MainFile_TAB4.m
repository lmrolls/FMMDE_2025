clear; clc; close;
%%#################################################################################
% CODE USED TO GENERATE TABLE 4
%
% NN_start: different cross section dimensions to loop over
% TT_start: different time series dimensions to loop over
% k0_start: different values of k0 to loop over
%
% rmax: max number of estimated factors,
% pval: critical value used in sequential testing procedure
% nreps: number of monte carlo replications
%
% mOutRAT, mOutTest: Output of Table 4, proportion of successful attempts
% in the estimation of the true number of factors rr in over the total
% number of monte carlo replications
%
% MAKES USE OF THE FUNCTION fLeeShaoNonLinear() which generates data from
% the factor model in EQ.5
%  THE FUNCTION GENERATES DATA FROM THREE LATENT FACTORS
%%#########################################################################################

addpath SW_MDDM
addpath numberF
addpath utils
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN_start    = [50 100 200];
TT_start    = [50 100 200];
k0_start = [ 1 10 25];

rmax  = 10;
pval  = 0.05;
nreps = 1000;

mOutRAT  = zeros(length(TT_start),length(NN_start),length(k0_start)); 
mOutTest = zeros(length(TT_start),length(NN_start),length(k0_start));

for k00 = 1:length(k0_start)
    k0 = k0_start(k00);
    for TT = 2:length(TT_start)
        T = TT_start(TT);
        disp(T)
        for NN = 1:length(NN_start)
            N = NN_start(NN);
            disp(N)
            
            mNfactors = zeros(nreps,2);
            parfor rep = 1:nreps
                disp(rep)
                x                = fLeeShaoNonLinear(T,N);
                mX               = standardize(x);    
                
                [~,~,~,~,icstar] = factorMDDM4(mX, k0, rmax);
                vPvals           = seqTest(mX,300,pval,k0);
                vNfactorsK0      = sum(vPvals<=pval);
                vNfactorsRatio   = icstar;
                
                
                
                mNfactors(rep,:) = [vNfactorsK0,vNfactorsRatio];
                
                
                
            end
            
            out      = mean(mNfactors==3,1);
            
            mOutTest(TT,NN,k00) = out(1);
            mOutRAT(TT,NN,k00) = out(2);
            
        end
    end
end


%%#########################################################################################
% Monte Carlo Simulation for Factor Model Estimation (Table 4)
%
% Purpose: 
%   Generates data from a nonlinear factor model (EQ.3 in the paper) with three latent factors
%   and evaluates the performance of two factor estimation methods: Eigenvalue 
%   Ratio and Sequential Testing.
%
% Parameters:
%   - NN_start: Cross-sectional dimensions [50, 100, 200]
%   - TT_start: Time-series dimensions [50, 100, 200]
%   - k0_start: Tuning parameter values [1, 10, 25]
%   - rmax: Maximum number of estimated factors (10)
%   - pval: Critical value for sequential testing (0.05)
%   - nreps: Number of Monte Carlo replications (100)
%
% Outputs:
%   - mOutRAT: Proportion of correct factor number estimates (3) using Eigenvalue Ratio
%   - mOutTest: Proportion of correct factor number estimates (3) using Sequential Testing
%
% Dependencies:
%   - fLeeShaoNonLinear(): Generates data from the factor model in EQ.3 in
%   the paper
%   - factorMDDM4(): Computes Eigenvalue Ratio-based factor estimation
%   - seqTest(): Performs sequential testing for factor estimation
%   - standardize(): Standardizes input data
%
% Notes:
%   - Random number generator set with rng(1, 'twister') for replicability
%   - Results stored in mOutRAT and mOutTest for Table 4
%
%%#########################################################################################

close all; clear all; clc;

addpath('../factorEstimation');
addpath('../DGPs');
addpath('subfunctions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set random number generator for replicability
rng(1, 'twister');


NN_start    = [50 100 200];
TT_start    = [50 100 200];
k0_start = [ 1 10 25];

rmax  = 10;
pval  = 0.05;
nreps = 300;

mOutRAT  = zeros(length(TT_start),length(NN_start),length(k0_start)); 
mOutTest = zeros(length(TT_start),length(NN_start),length(k0_start));

for k00 = 1:length(k0_start)
    k0 = k0_start(k00);
    for TT = 1:length(TT_start)
        T = TT_start(TT);
        disp(T)
        for NN = 1:length(NN_start)
            N = NN_start(NN);
            disp(N)
            
            mNfactors = zeros(nreps,2);
            parfor rep = 1:nreps
                disp(rep)

                substream = RandStream('mt19937ar', 'Seed', rep);
                RandStream.setGlobalStream(substream);

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

%%

%==================== =======  %
%                              %
%  EXCEL TABLES PRODUCTION     %
%                              %
%==============================%
filename = 'Tab4.xlsx';
sheet = 'Results';

% Headers
headers = {'', '', '', '', 'Sequential Test', '', '', '', 'Eigenvalue Ratio', '', ''};
subheaders = {'k0', '', 'T \ n', '', '50', '100', '200', '', '50', '100', '200'};

% Write to Excel
writecell(headers, filename, 'Sheet', sheet, 'Range', 'A3');
writecell(subheaders, filename, 'Sheet', sheet, 'Range', 'A4');

k0_values = k0_start;
T_values = TT_start;
row_start = 5;
for i = 1:length(k0_values)
    k00 = i; % Index for k0
    k0 = k0_values(i);
    for j = 1:length(T_values)
        TT = j; % Index for T
        T = T_values(j);
        
        % Extract data for this k0 and T, across all N
        seq_data = [mOutTest(TT, 1, k00), mOutTest(TT, 2, k00), mOutTest(TT, 3, k00)]; % N=50, 100, 200
        eigen_data = [mOutRAT(TT, 1, k00), mOutRAT(TT, 2, k00), mOutRAT(TT, 3, k00)]; % N=50, 100, 200
        
        % Round to 2 decimals
        seq_data = round(seq_data, 2);
        eigen_data = round(eigen_data, 2);
        
        % Compute row position
        row = row_start + (i-1)*4 + (j-1);
        
        % Write to Excel
        if j == 1
            % Write k0 only once per group, will merge later
            writecell({k0}, filename, 'Sheet', sheet, 'Range', sprintf('A%d', row));
        end
        writecell({T}, filename, 'Sheet', sheet, 'Range', sprintf('C%d', row));
        writematrix(seq_data, filename, 'Sheet', sheet, 'Range', sprintf('E%d', row)); % E:G
        writematrix(eigen_data, filename, 'Sheet', sheet, 'Range', sprintf('I%d', row)); % I:K
    end
    row_start = row_start + 1; % Add empty row between k0 values
end
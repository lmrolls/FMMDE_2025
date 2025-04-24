%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONTE CARLO SIMULATION FOR FACTOR MODEL ESTIMATION
%
% Purpose:
%   This script generates Tables 1 to 3 in the paper by simulating factor models 
%   under varying parameters (true number of factors, DGP, signal-to-noise ratio) 
%   and evaluating factor estimation accuracy using Sequential Testing and 
%   Eigenvalue Ratio methods.
%
% Key Functionality:
%   - Conducts Monte Carlo simulations with nested loops over:
%     1. True number of factors (r)
%     2. Data Generating Process (DGP)
%     3. Signal-to-noise parameter (theta)
%     4. Monte Carlo replications (via SimulFun)
%   - Stores results in tables for export, comparing estimation success rates.
%
% Parameters:
%   - n: Cross-sectional dimension (default: 100)
%   - T: Time-series dimension (default: 100)
%   - k0: Parameter for MDDM (default: 1)
%   - pval: Critical value for sequential testing (default: 0.05)
%   - nIters: Number of Monte Carlo replications (default: 200)
%   - r: True number of factors (2 to 5)
%   - DGP: Data Generating Process (1 to 3)
%   - thetaMethod: Signal-to-noise parameter multiplier
%       (1: 0.5*r, 2: r, 3: 3*r, 4: 5*r)
%
% Outputs:
%   - finalTable: Table with estimation success rates for Sequential Test and 
%                 Eigenvalue Ratio methods across parameter combinations.
%
% Notes:
%   - Requires 'factorEstimation' and 'subfunctions' folders in the MATLAB path.
%   - Random number generator seeded for replicability (rng(1, 'twister')).
%   - Parallelized Monte Carlo loop (parfor) for efficiency.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



close all; clear all; clc;

addpath('../factorEstimation');
addpath('subfunctions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set random number generator for replicability
rng(1, 'twister');

n      = 100;
T      = 100;
k0     = 1; 
pval   = 0.05; 

nIters = 200;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          EXPORTABLE TABLE CREATION COMMANDS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rvar=zeros(12,1);
dgp=zeros(12,1);
r=zeros(12,1);
rx5=zeros(12,1);
rx3=zeros(12,1);
rx05=zeros(12,1);

sequentialTest=table(rvar,dgp,rx05,r,rx3,rx5);
eigenvalueRatio=table(rvar,dgp,rx05,r,rx3,rx5);
finalTable=table(sequentialTest,eigenvalueRatio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           MONTE CARLO SIMULATION CODE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r   = 2;


while r <= 5           % number of factors that are simulated
    DGP = 1;
    while DGP <= 3     % DGP from 1 to 3 
        thetaMethod = 1;
        while thetaMethod <= 4   % Theta, the signal to noise parameter

            out = NaN(nIters,2);
            
            parfor k = 1:nIters  % number of Monte Carlo replications
                substream = RandStream('mt19937ar', 'Seed', k);
                RandStream.setGlobalStream(substream);
                
                out(k,:) = SimulFun(T,n,k0,r,thetaMethod,DGP,pval);
                disp([k,thetaMethod,DGP,r]);
                
            end
            disp(mean(out==r))
           

            index = (r-2)*3 + DGP; 
            
            finalTable.sequentialTest.rvar(index) = r;
            finalTable.sequentialTest.dgp(index)  = DGP;
            
    
            finalTable.eigenvalueRatio.rvar(index) = r;
            finalTable.eigenvalueRatio.dgp(index)  = DGP;
            


            if thetaMethod == 1

                finalTable.sequentialTest.rx05(index) = mean(out(:,1)==r);
                finalTable.eigenvalueRatio.rx05(index) = mean(out(:,2)==r);

            elseif thetaMethod ==2

                finalTable.sequentialTest.r(index) = mean(out(:,1)==r);
                finalTable.eigenvalueRatio.r(index) = mean(out(:,2)==r);

            elseif thetaMethod ==3

                finalTable.sequentialTest.rx3(index) = mean(out(:,1)==r);
                finalTable.eigenvalueRatio.rx3(index) = mean(out(:,2)==r);

            elseif thetaMethod ==4

                finalTable.sequentialTest.rx5(index) = mean(out(:,1)==r);
                finalTable.eigenvalueRatio.rx5(index) = mean(out(:,2)==r);

            end


            thetaMethod = thetaMethod + 1;
        end
       DGP = DGP + 1; 
    end
  r = r + 1;  
end

%%

%==================== =======  %
%                              %
%  EXCEL TABLES PRODUCTION     %
%                              %
%==============================%
filename = 'Simulation_Results.xlsx';
sheet = 'Results';

% Headers
headers = {'', '','','','', 'Sequential Test', '', '', '', '','', 'Eigenvalue Ratio', '', '', ''};
subheaders = {'r','',  'theta','','0.5r', 'r', '3r', '5r', '','', '0.5r', 'r', '3r', '5r'};

% Write to Excel
writecell(headers, filename, 'Sheet', sheet, 'Range', 'A3');
writecell(subheaders, filename, 'Sheet', sheet, 'Range', 'A4');

r_values = [2, 3, 4, 5];
dgp_labels = {'DGP1', 'DGP2', 'DGP3'};
row_start = 5;
for i = 1:length(r_values)
    r = r_values(i);
    for j = 1:length(dgp_labels)
        dgp = dgp_labels{j};
        row = row_start + (i-1)*4 + (j-1);
        index = (r-2)*3 + j; % Match simulation indexing
        
        % Extract data for this r and DGP
        seq_data = [finalTable.sequentialTest.rx05(index), ...
                    finalTable.sequentialTest.r(index), ...
                    finalTable.sequentialTest.rx3(index), ...
                    finalTable.sequentialTest.rx5(index)];
        eigen_data = [finalTable.eigenvalueRatio.rx05(index), ...
                      finalTable.eigenvalueRatio.r(index), ...
                      finalTable.eigenvalueRatio.rx3(index), ...
                      finalTable.eigenvalueRatio.rx5(index)];
        
        % Round to 2 decimals
        seq_data = round(seq_data, 2); % 0.5r, r, 3r, 5r for Sequential Test (C:G)
        eigen_data = round(eigen_data, 2); % 0.5r, r, 3r, 5r for Eigenvalue Ratio (I:L)
        
        % Write to Excel
        writecell({r}, filename, 'Sheet', sheet, 'Range', sprintf('A%d', row));
        writecell({dgp}, filename, 'Sheet', sheet, 'Range', sprintf('B%d', row)); % Theta column (DGP label)
        writematrix(seq_data, filename, 'Sheet', sheet, 'Range', sprintf('E%d', row)); % C:G (0.5r, r, 3r, 5r)
        writematrix(eigen_data, filename, 'Sheet', sheet, 'Range', sprintf('K%d', row)); % I:L (0.5r, r, 3r, 5r)
    end
    row_start = row_start + 1; % Add empty row between r values
end


%%#############################################################
% THIS FUNCTION IS USED TO 1) SIMULATE HIGH DIMENSIONAL DATA FROM
% PRESPECIFIED FACTOR MODEL DGPs 
% 2) ESTIMATE THE FACTOR MODELS AND NUMBER OF FACTORS
%
% MAKES USE OF THE FUNCTION BarigozziSim() TO SIMULATE DATA
% INPUTS OF THIS FUNCTION ARE
% 1,2) T,N: time and cross sectional dimensions
% 3)   r: true number of factors to be simulated
% 4)   theta: signal to noise ratio parameter type, from 1 to 4
% 5)   DGP: DGP type, from 1 to 3 as indicated in the paper
% 
% IN THE CASE OF THETA WE HAVE THE FOLLOWING
% 1 ==> theta = 0.5*r
% 2 ==> theta = r
% 3 ==> theta = 3*r
% 4 ==> theta = 5*r
%
% rmax is the maximum number of factors to be estimated
%%#######################################################

function[mNfactors,vPvals] = SimulFun(T,n,k0,r,theta,DGP,pval) 

mX = BarigozziSim(T,n,r,theta,DGP);
mX = standardize(mX);
rmax = 10;


[~,~,~,~,icstar] = factorMDDM(mX, k0, rmax);
vPvals=seqTest(mX,300,pval,k0); %  
vNfactorsK0      = sum(vPvals<=pval);
vNfactorsRatio   = icstar;


mNfactors = [vNfactorsK0,vNfactorsRatio];

end


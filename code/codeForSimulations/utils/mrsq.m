function [R2,mR2,t10_s,t10_mR2] = mrsq(Fhat,lamhat,series)
% =========================================================================
% DESCRIPTION
% This function computes the R-squared and marginal R-squared from
% estimated factors and factor loadings.
%
% -------------------------------------------------------------------------
% INPUTS
%           Fhat    = estimated factors (one factor per column)
%           lamhat  = factor loadings (one factor per column)
%           ve2     = eigenvalues of covariance matrix
%           series  = series names
%
% OUTPUTS
%           R2      = R-squared for each series for each factor
%           mR2     = marginal R-squared for each series for each factor
%           mR2_F   = marginal R-squared for each factor
%           R2_T    = total variation explained by all factors
%           t10_s   = top 10 series that load most heavily on each factor
%           t10_mR2 = marginal R-squared corresponding to top 10 series
%                     that load most heavily on each factor 
%           
% -------------------------------------------------------------------------
% NOTES
% Authors: Michael W. McCracken and Serena Ng
% Date: 6/7/2017
% Version: MATLAB 2014a
% Required Toolboxes: None
%
% =========================================================================
% FUNCTION

% N = number of series, ic = number of factors
[N,ic] = size(lamhat); 

% Preallocate memory for output 
R2 = NaN(N,ic);                           
mR2 = NaN(N,ic);
t10_s=cell(15,ic);
t10_mR2=NaN(15,ic);

% Compute R-squared and marginal R-squared for each series for each factor
for i = 1:ic
    R2(:,i)  = (var(Fhat(:,1:i)*lamhat(:,1:i)'))';  
    mR2(:,i) = (var(Fhat(:,i)*lamhat(:,i)'))';
end

% Sort series by marginal R-squared for each factor
[vals,ind] = sort(mR2,'descend');

% Get top 10 series that load most heavily on each factor and the
% corresponding marginal R-squared values
for i=1:ic
    t10_s(:,i)=series(ind(1:15,i));
    t10_mR2(:,i)=vals(1:15,i);
end


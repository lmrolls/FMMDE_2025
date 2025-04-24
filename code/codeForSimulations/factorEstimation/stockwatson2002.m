function [F,ehat,lambda,chat] = stockwatson2002(x,nfac)

% X is observed
% nfac is the estimated number of factors
% F is T by r matrix of true factors
% Lambda N by r is the true loading matrix
% C=F*Lambda' T by N is the true common component
% chat is the estimated common component


[~,bign]              =  size(x);
xx                    =  x'*x;
[Fhat0,eigval,Fhat1]  =  svd(xx);
lambda                =  Fhat0(:,1:nfac)*sqrt(bign);
F                     =  x*lambda/bign;
ehat                  =  x-F*lambda';
chat                  =  F*lambda';
ve2                   =  sum(ehat'.*ehat')'/bign;
ss                    =  diag(eigval);

end
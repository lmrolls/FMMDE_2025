
function[vPvals,out] = seqTestProva(Y,B,crit,k0,esta)   

if nargin == 4
  esta = "te";
end

r=15;

[F,Ahat,~] = factorMDDM(Y,k0,r);

[T,colsF] = size(F); vPvals = NaN(colsF,1);
 %colsF= floor(colsF/3);
 F = F(:,1:colsF); Ahat = Ahat(:,1:colsF);

if esta == "te"
  test  = @(x) testShao(x, 'w', 9, 0); 
  bootF = @(x) bernRadem(x);
elseif esta == "no"
  test  = @(x) testShao(x, 'w', 9, 0);
  bootF = @(x) bernEsc(x);
end

for i = 1:r   


        stat = test(F(:,i));  
        Tstar = zeros(B,1);

        for b=1:B
            W = bootF(T); Estar  = F(:,i:end);  

            for ik = 1:T
                Estar(ik,:) = Estar(ik,:) * W(ik);
            end

            Ystar       = [F(:,1:(i-1)),Estar]*Ahat';  Fstar = factorMDDM(Ystar,k0,r);
            Tstar(b)    = test(Fstar(:,i)); 

        end

        vPvals(i) = sum(Tstar>=stat )/(B+1);


        if vPvals(i) > crit
            for q = (i+1):colsF
                vPvals(q) = 1;
            end
            break            
        end
      



end

logik = vPvals <= crit;
out = F(:,logik);


end


%%


function[fhat,Ahat,chat] = factorMDDM(Y,k0,~)

[n, p] = size(Y);
S = zeros(p,p);                   
                                  
for k=1:k0                                    % k0:lags cumul. MDDM(pag 219)
    S = S + MDDM( Y(1:(n-k),:), Y((1+k):n,:) );       %cumulative MDDM(k0)
end

[eVec,~] = eig(S);                      %eigendec. cumulative MDDM(k0)

Ahat = eVec(:,(1:p));   


% Components
    fhat=Y*Ahat;
    
    chat=fhat*Ahat';
    
end



%%


function out = bernEsc(n)
  p = ( 1 + sqrt(5) )/( 2*sqrt(5));
  out = zeros(n,1);
  for i = 1:length(out)
    out(i) = binornd(1,p);
    out(i) = ( out(i) >= 1 )*(0.5)*(1-sqrt(5)) + ( out(i) < 1 )*(0.5)*(1+sqrt(5));
  end
  
end

%%

function [stat,pval] = testShao(Y, type, M, B)

if nargin == 1
    type = 'w';
    B = 999;
    M = 6;
elseif nargin == 2
    B = 999;
    M = 6;
elseif nargin == 3
    B = 999;
end

%norm(X,'fro')
[N, p] = size(Y);


if type == 'w'
    T = 0;
    for j=1:M
        omegaj = (N - j + 1)/(N * (j)^2);
        T = T + omegaj * norm( MDDM( Y(1:(N-j),:), Y((1+j):N,:) ),'fro' );
    end
    stat = N * T;
    
    
    if B ~= 0
        Tstar = zeros(B,1);
        parfor b = 1:B
            T = 0;
            Ystar = Y; W = bernRadem(N);
            
            for ik = 1:N
                Ystar(ik,:) = Ystar(ik,:)*W(ik);
            end
            
            for j=1:M
                omegaj = (N - j + 1)/(N * (j)^2);
                T = T + omegaj * norm( MDDM( Ystar(1:(N-j),:), Ystar((1+j):N,:) ),'fro' );
            end
            
            Tstar(b) = N * T;
            
            
        end
        pval = sum(Tstar>=stat )/(B+1);

        
    end

end



end
%%

function out = bernRadem(n)
  p = 0.5;
  out = zeros(n,1);
  for i = 1:length(out)
    out(i) = binornd(1,p);
    out(i) = ( out(i) >= 1 ) + ( out(i) < 1 )*(-1);
  end
  
end

%%

function [mMDDM] = MDDM(mX, mY)  

% Computes the Martingale Difference Divergence Matrix (MDDM) between two
% vector time series X and Y, as defined in:
% Lee, C. E., & Shao, X. (2018). Martingale Difference Divergence Matrix and
% Its Application to Dimension Reduction for Stationary Multivariate Time Series.

 [cn, cp] = size(mX);
[~, cq] = size(mY);
mMDDM = zeros(cq,cq);
cYbar = mean(mY, 1);
mY = mY - ones(cn,1)*(cYbar);


if cp == 1
    for i=1:cn
        cXdist = abs(mX(i) - mX);
        mMDDM = mMDDM + mY(i,:) * (transpose(cXdist) * mY);
    end
    
else
    for i=1:cn
        
        cXdist = transpose(transpose(mX(i,:)) - transpose(mX));
        cXdist = cXdist.^2;
        cXdist = sqrt(sum(cXdist,2));
        mMDDM = mMDDM + transpose(mY(i,:)) * (transpose(cXdist) * mY);
    end
    
    
end
    mMDDM = -mMDDM/(cn^2);
    
end








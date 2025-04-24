
function[vPvals,out] = seqTest(Y,B,crit,k0,esta)   

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
  test  = @(x) escancianoMart2(x,0);
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


function[fhat,Ahat,chat] = factorMDDM(Y,k0,r)

[n, p] = size(Y);
S = zeros(p,p);                   
                                  
for k=1:k0                                    % k0:lags cumul. MDDM(pag 219)
    S = S + MDDM( Y(1:(n-k),:), Y((1+k):n,:) );       %cumulative MDDM(k0)
end

[eVec,eVal] = eig(S);                      %eigendec. cumulative MDDM(k0)
% eVal = diag(eVal);
% eVal = sort(eVal,'descend');
% 
% if p < 100
%     R  = p-1;
% elseif p>=100 && p<400
%     R = (p/2);
% else
%     R = floor((p/3));
% end
% 
% lambda = zeros(R,1);
% 
% for i=1:R
%     lambda(i) = eVal(i+1)/eVal(i);         %lambda= ratios subsequent eigenVals
% end
% 
% 
% [~,argMin] =  min(lambda(1:R));                 % argmin of ratios 
%    icstar = (argMin);   

% [~,argMin] =  min(lambda(2:R));                 % argmin of ratios 
%    icstar = (argMin+1);                      %for estimating size factor process                              
Ahat = eVec(:,(1:p));   


% Components
    fhat=Y*Ahat;
    
    chat=fhat*Ahat';
    %ic1 = argMin;
    %ss = eVal;
    %icstar = (argMin+2);
end


%%

function[stat,pval] = escancianoMart2(x, B)
[n] = length(x);
sigma2 = var(x);
xbar = mean(x);
stat3 = 0;
alfii = cell(n-1,1);

for j = 1:(n-1)
    nj = (n - j + 1);
    coeff = 1/(sigma2*((j*pi)^2)*nj);
    stat2 = 0;
    alfio = zeros(length((j+1) : n),length((j+1) : n));
    for t = (j+1) : n
        in  = x(t)-xbar;
        in2 = x(t-j);
        stat1 = 0;
        for s = (j+1) : n
            alfio2 = exp(-0.5*(in2-x(s-j))^2);
            stat1 = stat1 + in * ( x(s) - xbar ) * alfio2;
            alfio(t-j,s-j) = alfio2; 
            
            
        end
        stat2 = stat2 + stat1;
    end
    stat3 = stat3 + coeff * stat2;
    alfii{j} = alfio;
    
end

stat = stat3;

if B ~= 0
    Tstar = zeros(B,1);
    for b = 1:B
        stat3 = 0;
        Ystar = x; W = bernEsc(n);
        for ik = 1:n
            Ystar(ik,:) = Ystar(ik,:)*W(ik);
        end
        Ybar = mean(Ystar);
        sigma2star = var(Ystar);
        
        for j = 1:(n-1)
            alfio = alfii{j};
            nj = (n - j + 1);
            coeff = 1/(sigma2star*((j*pi)^2)*nj);
            stat2 = 0;

            for t = (j+1) : n
                in  = Ystar(t)-Ybar;
                %in2 = x(t-j);
                stat1 = 0;
                for s = (j+1) : n
                    stat1 = stat1 + in * (Ystar(s)-Ybar) * alfio(t-j,s-j);
                    
                end
                stat2 = stat2 + stat1;
            end
            stat3 = stat3 + coeff * stat2;
            
        end
        
        Tstar(b) = stat3;
        
        
    end
    
    pval = sum(Tstar>=stat )/(B+1);
   
    
else
    pval=2;
end


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



function[ip] = d_inner(X,Y)

[n,~] = size(X);
ip = sum(X .* Y,'all') ./ n ./ n ;


end

%%

function[DA] = d_center(A)

[n,~]=size(A);
R = sum(A,2);
C = sum(A,1);
S = sum(A,"all");
r = repmat(R,1,n)./(n);
c = transpose(repmat(C',1,n))./(n);
t = (zeros(n,n)+S)./(n^2);

 DA = A - r - c + t;
 %UA = UA - diag(diag(UA));
  
 clear A
 clear r
 clear c
 clear t
 
end



%%


function [Dx,Dy] = euclideanDist(mX,mY)

% Distance covariance. X is n times K


if nargin == 1
   mY = mX; 
end

[cn , ~] = size(mX); 
mD2X = zeros(cn,cn); 
mQX = mX * mX'; 

mD2Y = zeros(cn,cn); 
mQY = mY * mY'; 

for i=1:cn
    for j=(i+1):cn
        
        mD2X(i,j) = mQX(i,i) + mQX(j,j) - 2 * mQX(i,j); mD2X(j,i) = mD2X(i,j);
        mD2Y(i,j) = mQY(i,i) + mQY(j,j) - 2 * mQY(i,j); mD2Y(j,i) = mD2Y(i,j);
        
    end
end
% mDX = sqrt(mD2X); mDY = sqrt(mD2Y); 
Dx = sqrt(mD2X);
Dy = sqrt(mD2Y);
clear mD2X
clear mD2Y
clear mQX
clear mQY

end


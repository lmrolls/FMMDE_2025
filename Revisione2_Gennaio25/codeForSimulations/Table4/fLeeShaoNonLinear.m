function Yt = fLeeShaoNonLinear(T,N)

if nargin == 2
    r    = 3;
    A          = zeros(N,3);
    A(1:round(N/2),:) = -2 + (2+2)*rand(round(N/2),r);
    
end

burn    = 500;
T       = T+burn;
zt      = randn(T,1);
wt(1)   = 0;
for i = 2:T
    if wt(i-1,1) < 5
        wt(i,1) = 0.5 + (0.05*exp(-0.01*wt(i-1,1).^2) +0.9)*wt(i-1,1) + zt(i,:);
    else
        wt(i,1) = (0.9*exp(-10*wt(i-1,1).^2))*wt(i-1,1) + zt(i,:); 
    end
end

xt      = lagmatrix(wt,1:(r-1));
xt      = [wt,xt];
xt      = xt(burn+1:end,:);
sigma   = zeros(N,N);
H       = 0.9;
for i = 1:N
    for j = 1:N
        ind = i-j;
        if ind == 0
            sigma(i,j) = 1;
        else
        sigma(i,j) = 0.5*((abs(i-j)+1).^(2*H)-2*abs(i-j).^(2*H)+(abs(i-j)-1).^(2*H));
        end
    end
end
%et         = mvnrnd(zeros(N,1),eye(N),T);
et         = mvnrnd(zeros(N,1),0.25*sigma,T);
et         = et(burn+1:end,:);

Yt         = xt*A' + et;
function Yt = fLeeShao(T,N,r)
burn    = 500;
T       = T+burn;
%zt      = randn(T,1);
%b       = 1;
%a       = [1 -0.2];
%wt      = filter(b,a,zt);


Mdl = arima('Constant',0,'MA',{0.2},'Variance',1);
wt = simulate(Mdl,T);

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
et         = mvnrnd(zeros(N,1),sigma,T);
et         = et(burn+1:end,:);
A          = zeros(N,3);
A(1:round(N/2),:) = -2 + (2+2)*rand(round(N/2),r);
%A=  -2 + (2+2)*rand(N,r);
Yt         = xt*A' + et;
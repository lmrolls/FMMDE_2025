function[mX] = BarigozziSim(T,n,r,thetaMethod,DGP)

if thetaMethod == 1
    theta = 0.5 * r;
    
elseif thetaMethod == 2
    theta = r;
    
elseif thetaMethod == 3
    theta = 3 * r;
    
elseif thetaMethod == 4
    theta = 5 * r;

end



if DGP == 1
    
    %Homoskedastic idiosyncratic components
    mEpsi = randn(T,n);
    
elseif DGP == 2
    
    %Heteroskedastic idiosyncratic components
    mEpsi = zeros(T,n);
    for t = 1:T
        mEpsi(t,:) = (mod(t,2)~=0)*(randn(1,n)) +  (mod(t,2)==0)*(randn(1,n) + randn(1,n));
    end
    
elseif DGP == 3
    
    %Cross-sectional correlations among idiosyncratic components
    
    J     = max(n/20,10); beta = 0.2;
    vV    = randn(T, n+2*J);
    mBeta = zeros(n+2*J, n);
    vBeta = [beta*ones(J,1);1;beta*ones(J,1)];
    span  = length(vBeta);
    
    for k = 1:n
        mBeta(k:(span+k-1),k) = vBeta;
    end
    
    mEpsi = vV * mBeta;
    
elseif DGP == 4
    
    
    mEpsi = zeros(T,n); for i = 1:n
       
        par1 = 2; par2 = 2;
        while(par1+par2 >=1)
            par1=unifrnd(0.0,0.9);
            par2=unifrnd(0.0,0.9);           
        end
        
        Mdl = garch('Constant',0.0001,'GARCH',par1,'ARCH',par2);
        [~,mEpsi(:,i)] = simulate(Mdl,T,'NumPaths',1);
        
    end
   
    
   
    
end

mLambda = randn(r,n);
mF = zeros(T,r);



for i = 1:r 
par=unifrnd(0.5,0.9);
 e=randn(T,1);
 b = 1;
 a = [1 -par];
 mF(:,i) = filter(b,a,e);

end



mX = mF * mLambda + sqrt(theta) * mEpsi;
end


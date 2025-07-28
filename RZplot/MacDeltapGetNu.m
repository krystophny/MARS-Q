% calculate nu and Delta from giving function y(x)
% using order-N expansion to best fit the full data

function [nu,Delta,Fmin] = MacDeltapGetNu(x1,y1,N1)

global N NN K mi JJ
global x y A  

N     = N1;
[K,M] = size(y1);
NN    = N+N;
x     = x1;
y     = y1;

if NN>K  
   disp('MacDeltapGetNu: not enough data points!')
end

nu    = zeros(M,1);
Delta = zeros(M,1);
Fmin  = zeros(M,1);

nu_min =  0.01;
nu_max =  0.99;

%J  = floor(K/NN);
%JJ = [0:NN-1]*J + 1;
JJ = 1:NN;

x  = x(:);
A  = zeros(NN,NN);

for mi=1:M
    [nu(mi),Fmin(mi)] = fminbnd(@MacDeltapFun,nu_min,nu_max);
end

Delta = MacDeltapPlot(x,y,nu,N,JJ,1);


  
        



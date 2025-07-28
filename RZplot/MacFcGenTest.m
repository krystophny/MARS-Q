% generate excitation matrix in terms of Fourier harmonics

M = 6;
t = 2*pi/M;

m = [1:M]';
Y = zeros(M,M);
for k=0:M-1
  Y(:,k+1) = exp(i*m*k*t);
end

CondY = cond(Y)

X = inv(Y);
CondX = cond(X)

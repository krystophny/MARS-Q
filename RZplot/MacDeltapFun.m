% the objective function for the optimiser 
% called from MacDeltapGetNu.m

function fun = MacDeltapFun(nui)

global N NN K mi JJ
global x y A 

xn = x(JJ);
b  = y(JJ,mi);
for ni=1:N
    A(:,ni) = xn.^(ni-nui);
    A(:,ni+N) = xn.^(ni+nui);
end
u  = inv(A)*b;
ye = y(:,mi);
yn = zeros(K,1);
for ni=1:N
    yn = yn+u(ni)*x.^(ni-nui);
    yn = yn+u(ni+N)*x.^(ni+nui);
end

fun = sqrt(sum((yn-ye).^2))/sqrt(sum(ye.^2));




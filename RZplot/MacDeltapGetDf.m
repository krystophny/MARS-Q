% get the large solution coefficients
% from the full solution (KDELTAP=6 from MacReadBPLASMA.m)
% using LS fitting 
% y = c_1 + c1_2*|x|^(-nui)

function Delta = MacDeltapGetDf(x,y,nui,JJ,kplot)

NN    = length(JJ);
[K,M] = size(y);
N     = 2;

if N>NN, disp('MacDeltapGetDf: too few data points!'); end

Delta = zeros(M,1);
xn    = x(JJ);
xnn   = abs(xn).^(-nui);
b     = zeros(N,1);
A     = zeros(N,N);

if kplot==1, hf=figure; end

for mi=1:M
    b(1)   = sum(y(JJ,mi));
    b(2)   = sum(y(JJ,mi).*xnn);
    A(1,1) = NN;
    A(1,2) = sum(xnn);
    A(2,1) = A(1,2);
    A(2,2) = sum(xnn.^2);

    u  = inv(A)*b;
    Delta(mi) = u(2);

    if kplot==1 & mi==4
       if mi==4, res_u = u, end
       yn = u(1); 
       yn = yn + u(2)*abs(x).^(-nui);
       plot(x,y(:,mi),'r+',x,yn,'b-','LineWidth',2); hold on,
    end
end

if kplot==1
   xlabel('s-s_m','FontSize',16,'FontWeight','Bold')
   ylabel('b_m^1','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
end
   

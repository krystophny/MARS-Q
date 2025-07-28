% get the small solution coefficients
% using LS fitting for two terms 

function Delta = MacDeltapGetDls(x,y,nui,JJ,kplot)

N     = length(JJ);
[K,M] = size(y);

Delta = zeros(M,1);
xn    = x(JJ);
A     = zeros(2,2);
b     = zeros(2,1);

if kplot==1, hf=figure; end

for mi=1:M
    b(1)   = sum(y(JJ,mi).*xn.^(-nui));
    b(2)   = sum(y(JJ,mi).*xn.^(nui));
    A(1,1) = sum(xn.^(-2*nui));
    A(1,2) = N;
    A(2,1) = N;
    A(2,2) = sum(xn.^(2*nui));
    u  = inv(A)*b;
    Delta(mi) = u(2);

    if kplot==1 
       yn = u(1)*x.^(-nui) + u(2)*x.^nui;
       plot(x,y(:,mi),'r+',x,yn,'b-','LineWidth',2); hold on,
    end
end

if kplot==1
   xlabel('s-s_m','FontSize',16,'FontWeight','Bold')
   ylabel('\partial b_m^1/\partial s','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
end
   

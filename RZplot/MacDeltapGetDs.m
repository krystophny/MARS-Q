% get the small solution coefficients
% using LS fitting 
% y = c_1*x^(1+nui) + c_2*x^(2+nui) + ...

function Delta = MacDeltapGetDs(x,y,nui,N,JJ,kplot)

NN    = length(JJ);
[K,M] = size(y);

if N>NN, disp('MacDeltapGetDs: too few data points!'); end

Delta = zeros(M,1);
xn    = x(JJ);
b     = zeros(N,1);
A     = zeros(N,N);

if kplot==1, hf=figure; end

for mi=1:M
    for ki = 1:N
        b(ki)   = sum(y(JJ,mi).*xn.^(ki+nui));
        for ji=1:N
            A(ki,ji) = sum(xn.^(ki+ji+2*nui));
        end
    end
    u  = inv(A)*b;
    Delta(mi) = u(1);
    
    if kplot==1 
       yn = 0;
       for ki=1:N
           yn = yn + u(ki)*x.^(ki+nui);
       end
       plot(x,y(:,mi),'r+',x,yn,'b-','LineWidth',2); hold on,
    end
end

if kplot==1
   xlabel('s-s_m','FontSize',16,'FontWeight','Bold')
   ylabel('b_m^1','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
end
   

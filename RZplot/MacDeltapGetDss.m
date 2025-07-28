% get the small solution coefficients
% using LS fitting 
% y = c_1*x + c_2*x^(1+nui) + c_3*x^(2+nui) + ...

function Delta = MacDeltapGetDss(x,y,nui,N,JJ,kplot)

NN    = length(JJ);
[K,M] = size(y);

if N+1>NN, disp('MacDeltapGetDss: too few data points!'); end

Delta = zeros(M,1);
xn    = x(JJ);
b     = zeros(N+1,1);
A     = zeros(N+1,N+1);

if kplot==1, hf=figure; end

for mi=1:M
    b(1)   = sum(y(JJ,mi).*xn);
    A(1,1) = sum(xn.^2);
    for ki = 1:N
        A(1,ki+1) = sum(xn.^(ki+nui+1));
        b(ki+1)   = sum(y(JJ,mi).*xn.^(ki+nui));
        A(ki+1,1) = A(1,ki+1);
        for ji=1:N
            A(ki+1,ji+1) = sum(xn.^(ki+ji+2*nui));
        end
    end
    u  = inv(A)*b;
    Delta(mi) = u(2);

    
    if mi==4
       res_u = u
    end

    if kplot==1 
       yn = u(1)*x;
       for ki=1:N
           yn = yn + u(ki+1)*x.^(ki+nui);
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
   

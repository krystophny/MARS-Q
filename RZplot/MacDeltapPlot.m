% plot the optimal solution

function Delta = MacDeltapPlot(x,y,nu,N,JJ,kplot)

NN = N+N;
[K,M] = size(y);

Delta = nu;
xn    = x(JJ);

if kplot==1, hf=figure; end

for mi=1:M
    nui = nu(mi);
    b   = y(JJ,mi);
    for ni=1:N
        A(:,ni) = xn.^(ni-nui);
        A(:,ni+N) = xn.^(ni+nui);
    end
    u  = inv(A)*b;
    Delta(mi) = u(N+1);

    if kplot==1 
       yn = zeros(K,1);
       for ni=1:N
           yn = yn+u(ni)*x.^(ni-nui);
           yn = yn+u(ni+N)*x.^(ni+nui);
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
   

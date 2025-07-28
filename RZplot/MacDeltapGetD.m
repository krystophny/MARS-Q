% get the small solution coefficients
% assuming a portion of large solution is present

function Delta = MacDeltapGetD(x,y,nui,N,JJ,kplot)

NN = N+N;
[K,M] = size(y);

Delta = zeros(M,1);
xn    = x(JJ);
A     = zeros(NN,NN);

if kplot==1, hf=figure(1); end

for mi=1:M
    b   = y(JJ,mi);
    for ni=1:N
        A(:,ni) = xn.^(ni-nui-1);
        A(:,ni+N) = xn.^(ni+nui-1);
    end
    u  = inv(A)*b
    Delta(mi) = u(N+1);

    if kplot==1 
       yn = zeros(K,1);
       for ni=1:N
           if ni>1, yn = yn+u(ni)*x.^(ni-nui-1); end
           yn = yn+u(ni+N)*x.^(ni+nui-1);
       end

       %u(1) = -0.25;
       yn = y(:,mi) - u(1)*x.^(-nui) - u(2)*x.^(nui);
       %plot(x,y(:,mi),'r+',[0;x],[0;yn],'b-+','LineWidth',2); hold on,
       %loglog(x,abs(yn),'b-+'), hold on
    end
end

if kplot==1
	plot(x,y,'-+','LineWidth',2); hold on,
   xlabel('s-s_m','FontSize',16,'FontWeight','Bold')
   ylabel('\partial b_m^1/\partial s','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
end
   

% get the large solution coefficients
% from the full solution (KDELTAP=6 from MacReadBPLASMA.m)
% using LS fitting 
% y = c_1*|x|^(-nu) + c_2*|x|^(1-nu) + c_3*|x|^(1+nu) + ...
%     c_4*x + c_5

function Delta = MacDeltapGetDff(x,y,nui,JJ,kplot)

NN    = length(JJ);
[K,M] = size(y);
N     = 5;

if N>NN, disp('MacDeltapGetDff: too few data points!'); end

Delta = zeros(M,1);
xn    = x(JJ);
b     = zeros(N,1);
A     = zeros(N,N);

if kplot==1, hf=figure; end

for mi=1:M
    b(1)   = sum(y(JJ,mi).*abs(xn).^(-nui));
    b(2)   = sum(y(JJ,mi).*abs(xn).^(1-nui));
    b(3)   = sum(y(JJ,mi).*abs(xn).^(1+nui));
    b(4)   = sum(y(JJ,mi).*xn);
    b(5)   = sum(y(JJ,mi));

    A(1,1) = sum(abs(xn).^(-2*nui));
    A(1,2) = sum(abs(xn).^(1-2*nui));
    A(1,3) = sum(abs(xn));
    A(1,4) = sum(abs(xn).^(-nui).*xn);
    A(1,5) = sum(abs(xn).^(-nui));

    A(2,1) = A(1,2);
    A(2,2) = sum(abs(xn).^(2-2*nui));
    A(2,3) = sum(abs(xn).^2);
    A(2,4) = sum(abs(xn).^(1-nui).*xn);
    A(2,5) = sum(abs(xn).^(1-nui));

    A(3,1) = A(1,3);
    A(3,2) = A(2,3);
    A(3,3) = sum(abs(xn).^(2+2*nui));
    A(3,4) = sum(abs(xn).^(1+nui).*xn);
    A(3,5) = sum(abs(xn).^(1+nui));

    A(4,1) = A(1,4);
    A(4,2) = A(2,4);
    A(4,3) = A(3,4);
    A(4,4) = sum(xn.^2);
    A(4,5) = sum(xn);

    A(5,1) = A(1,5);
    A(5,2) = A(2,5);
    A(5,3) = A(3,5);
    A(5,4) = A(4,5);
    A(5,5) = NN;

    u  = inv(A)*b;
    Delta(mi) = u(1);

    if kplot==1 & mi==4
       if mi==4, res_A=A, res_u=u, end
       yn =      u(1)*abs(x).^(-nui); 
       yn = yn + u(2)*abs(x).^(1-nui);
       yn = yn + u(3)*abs(x).^(1+nui);
       yn = yn + u(4)*x;
       yn = yn + u(5);
       plot(x,y(:,mi),'r+',x,yn,'b-','LineWidth',2); hold on,
    end
end

if kplot==1
   xlabel('s-s_m','FontSize',16,'FontWeight','Bold')
   ylabel('b_m^1','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
end
   

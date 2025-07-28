% Plot plasma rotation induced resonance splitting 
% with parallel Alfven and sound wave continua 
% Need to read the following data from MARS-F output
%   PROFEQ.OUT   
% Note that this is a stand-alone script w/o calling other Mac-procedures

% input parameter
GAM = 5/3;  %ratio of specfic heats
kfs = 0;    %=1: plot also fast and slow magneto-acoustic freqeuncies

% read in data
data = load([SDIR 'PROFEQ.OUT']);
s    = data(:,1);
q    = data(:,2);
P    = data(:,4);
rho  = data(:,5);
w    = data(:,6);

% find rational surfaces
[Is,rs,qqs] = MacGetRatSurf(s.^2,q,n,[0:100],9);
NS = length(Is);
rs = sqrt(rs);
mm = round(abs(n)*qqs);

% find splitted resonant surfaces
% due to both shear Alfven and parallel sound waves
% plot also fast and slow magneto-acoustic waves

hf=figure(kplot_cont);

% get q-intervals
qq=min(q);
for m=1:NS-1
    qq = [qq qqs(m)+0.5];
end
qq = [qq max(q)];

% plot all frequencies in each q-interval
for k=1:NS
    II = find(q>=qq(k) & q<=qq(k+1));
    s1 = s(II);
    q1 = q(II);
    P1 = P(II);
    w1 = w(II);
    rho1 = rho(II);
    m   = mm(k);

    if kaxisx==1, x1 = s1; xlab = '\psi_p^{1/2}'; end
    if kaxisx==2, x1 = s1.^2; xlab = '\psi_p'; end
    if kaxisx==3, x1=q1; xlab='q'; end

    Fps = 1 + (q1./(m-1-n*q1)).^2 + (q1./(m+1-n*q1)).^2;
    P2  = P1.*Fps;
    q11 = (m./q1-n);
    eps = s1/ASPECT;
    a2  = 4*q11.^2*GAM.*P2./(1+GAM*P2).^2./(n^2+m^2./eps.^2);
    wa2 = q11.^2./Fps./rho1;
    wh2 = GAM*P2./(1+GAM*P2).*wa2;
    wf2 = 0.5*(1+GAM*P2).*(n^2+m^2./eps.^2)./rho1./Fps.*(1+sqrt(1-a2));
    ws2 = 0.5*(1+GAM*P2).*(n^2+m^2./eps.^2)./rho1./Fps.*(1-sqrt(1-a2));

    semilogy(x1,abs(w1),'k-','LineWidth',2), hold on,
    semilogy(x1,sqrt(wa2),'r-',x1,sqrt(wh2),'b-','LineWidth',1), hold on,
    if (kfs==1) 
       semilogy(x1,sqrt(wf2),'r--',x1,sqrt(ws2),'b--','LineWidth',0.5), hold on,
    end

    xlabel(xlab,'FontSize',18,'FontWeight','Bold') 
    ylabel('frequencies normalised by \omega_A','FontSize',16,'FontWeight','Bold')
    ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
    if (kfs==1)
       legend('|\omega|','|\omega_a|','|\Omega_h|','|\omega_f|','|\Omega_s|') 
    else    
       legend('|\omega|','|\omega_a|','|\Omega_h|') 
    end
end

% calculate location of rational surfaces
% giving the q(s) and the toroidal mode number n
% this function works also for non-monotonic q-profile

function [Ires,Sres,Qres]=MacGetRatSurf(s,q,n,mm,kplot,SCOL)
 
qmin = min(q);
qmax = max(q);
mmin = max([ceil(abs(n)*qmin)  min(mm)]);
mmax = min([floor(abs(n)*qmax) max(mm)]);
mn   = [mmin:mmax]';

Ires = [];
Sres = [];
Qres = [];

if abs(n)>0
for k=1:length(s)-1
    II = find( (q(k)*abs(n)-mn).*(q(k+1)*abs(n)-mn)<=0 );
    if length(II)>0
       Ires = [Ires; k*ones(length(II),1)];
       Qres = [Qres; mn(II)/abs(n)];
    end
end
 
Sres = s(Ires);

end

if kplot>0
 
hf = figure(kplot);
plot(s.^2,q,'LineWidth',3,'Color',SCOL), hold on,
for k=1:length(Ires)
    plot([Sres(k)^2 Sres(k)^2],[0 Qres(k)],'k--')
    plot([0 Sres(k)^2],[Qres(k) Qres(k)],'k--')
end
xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
ylabel('q','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

end

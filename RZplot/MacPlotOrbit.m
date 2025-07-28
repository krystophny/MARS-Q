%plot particle orbit computed by MARS-K
%and saved to ORBIT.OUT

function MacPlotOrbit(filename)

global Mac

data = load(filename);

ORBS0 = data(1,1);
ORBQ0 = data(1,2);
ORBQP = data(1,3);

ORBL  = data(2,2);
ORBSIG= data(2,3);
ORBTB = data(2,4);
ORBWD = data(2,5);
ORBWTH= data(2,6);

t     = data(3:end,1);
s     = data(3:end,2);
chi   = data(3:end,3);
phi   = data(3:end,4); 
Rorb  = data(3:end,5);
Zorb  = data(3:end,6);

if Mac.plot_shape > 0
   figure(Mac.plot_shape)
   plot(Rorb*Mac.R0EXP,Zorb*Mac.R0EXP,'b-'), hold on,
end

hf = figure(Mac.plot_orbit*10 + 1);
plot(t,s,'b-','LineWidth',2), hold on,
plot([t(1) t(end)],[ORBS0 ORBS0],'k--')
plot([t(1) t(end)],[max(s) max(s)],'b--')
plot([t(1) t(end)],[max(s) max(s)]-ORBWTH,'b--')
xlabel('time/\tau_A','FontSize',18,'FontWeight','Bold'),
ylabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold'),
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),

hf = figure(Mac.plot_orbit*10 + 2);
plot(t,chi,'b-','LineWidth',2), hold on,
xlabel('time/\tau_A','FontSize',18,'FontWeight','Bold'),
ylabel('\chi','FontSize',18,'FontWeight','Bold'),
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),

hf = figure(Mac.plot_orbit*10 + 3);
plot(t,phi,'r-','LineWidth',2), hold on,
SL{1} = 'total';
if ORBL>0.5  %trapped particle
   SL{2} = 'bounce';
   plot(t,t*0,'b-','LineWidth',2)
   SL{3} = 'precession';
   plot(t,t*ORBWD,'b--','LineWidth',2)
   phin = phi-t*ORBWD;
else         %passing
   SL{2} = 'transit';
   plot(t,t*2*pi/ORBTB*ORBQ0*ORBSIG,'b-','LineWidth',2)
   SL{3} = 'precession';
   plot(t,t*ORBWD,'b--','LineWidth',2)
   phin = phi-t*2*pi/ORBTB*ORBQ0*ORBSIG - t*ORBWD*1;
end
SL{4} = 'residual';   
plot(t,phin,'k-','LineWidth',2), hold on,
plot([t(1) t(end)],[max(phin) max(phin)],'k--'),
plot([t(1) t(end)],[min(phin) min(phin)],'k--'),
xlabel('time/\tau_A','FontSize',18,'FontWeight','Bold'),
ylabel('\phi','FontSize',18,'FontWeight','Bold'),
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),
legend(SL)



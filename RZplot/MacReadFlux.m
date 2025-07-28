% read and compare MARS-Q computed particle fluxes
% due to MHD term and NTV

kdev    = 25; 

MacReadParam

SDIR = '~/IRIS/WorkD3Da/';

a = load([SDIR 'TORQUEREY.OUT']);
b = load([SDIR 'TORQUENTV.OUT']);
c = load([SDIR 'PROFEQ.OUT']);

%s     = (c(1:end-1,2)+c(2:end,2))/2;
s     = b(:,1);
Fntv  = b(:,3);
Fntvi = b(:,4);
Fntve = b(:,5);
Fmhd  = a(:,3);

hf = figure(1);
plot(s,Fmhd,'r-','LineWidth',3), hold on,
plot(s,Fntv,'b-','LineWidth',3), hold on,
plot(s,Fntvi,'b--','LineWidth',2), hold on,
plot(s,Fntve,'b-.','LineWidth',2), hold on,
xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
%xlabel('q','FontSize',18,'FontWeight','Bold')
ylabel('Particle Flux','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
legend('MHD','NTV','NTVi','NTVe')
%a=axis;
%plot([11/3 11/3],[a(3) a(4)],'k-')


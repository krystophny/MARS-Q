% read data METRICS.OUT and plot 

SDIR  = '/.automount/funsrv4.ccfe.ac.uk/root/home1/yliu/MAST_ADD/';
plotB = 2;
 
% reading data
data   = load([SDIR 'METRICS.OUT']);
N      = size(data,1);
II     = [2:N];
s      = data(II,1);
DG11m0 = data(II,2)+data(II,3)*i;
DG12m0 = data(II,4)+data(II,5)*i;
DG22m0 = data(II,6)+data(II,7)*i;
DG33m0 = data(II,8)+data(II,9)*i;
DG11m1 = data(II,10)+data(II,11)*i;
DG12m1 = data(II,12)+data(II,13)*i;
DG22m1 = data(II,14)+data(II,15)*i;
DG33m1 = data(II,16)+data(II,17)*i;

% plottling
hf=figure(10*plotB + 1);
plot(s,real(DG11m0),'r-',s,imag(DG11m0),'r--','LineWidth',2), hold on,
plot(s,real(DG11m1),'b-',s,imag(DG11m1),'b--','LineWidth',2), hold on,
xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
ylabel('g_{11}/J','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('Re(m=0)','Im(m=0)','Re(m=1)','Im(m=1)')

hf=figure(10*plotB + 2);
plot(s,real(DG12m0),'r-',s,imag(DG12m0),'r--','LineWidth',2), hold on,
plot(s,real(DG12m1),'b-',s,imag(DG12m1),'b--','LineWidth',2), hold on,
xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
ylabel('g_{12}/J','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('Re(m=0)','Im(m=0)','Re(m=1)','Im(m=1)')

hf=figure(10*plotB + 3);
plot(s,real(DG22m0),'r-',s,imag(DG22m0),'r--','LineWidth',2), hold on,
plot(s,real(DG22m1),'b-',s,imag(DG22m1),'b--','LineWidth',2), hold on,
xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
ylabel('g_{22}/J','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('Re(m=0)','Im(m=0)','Re(m=1)','Im(m=1)')

hf=figure(10*plotB + 4);
plot(s,real(DG33m0),'r-',s,imag(DG33m0),'r--','LineWidth',2), hold on,
plot(s,real(DG33m1),'b-',s,imag(DG33m1),'b--','LineWidth',2), hold on,
xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
ylabel('g_{33}/J','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('Re(m=0)','Im(m=0)','Re(m=1)','Im(m=1)')


% plotting TMP_METRIC
if 1==0
data   = load([SDIR 'TMP_METRIC.OUT']);
NP     = data(1,1);
M      = data(1,2);
NV     = data(NP*M+2,1);
N      = NP + NV;
II     = [2:NP*M+1 NP*M+3:NP*M+NV*M+2];
G11    = data(II,1);
G22    = data(II,2);
JAC    = data(II,3);

G11    = reshape(G11,M,N);
G22    = reshape(G22,M,N);
JAC    = reshape(JAC,M,N);

G11 = G11(:,2:end);
G22 = G22(:,2:end);
JAC = JAC(:,2:end);

hf=figure(10*plotB + 5);
plot(s,sum(G11,1),'b-',s,sum(G22,1),'k-',s,sum(JAC,1),'r-','LineWidth',2), hold on,
xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
ylabel('y','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('g_{11}','g_{22}','J')

hf=figure(10*plotB + 6);
plot(s,G11), hold on,
xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
ylabel('G11','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf=figure(10*plotB + 7);
plot(s,G22), hold on,
xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
ylabel('G22','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf=figure(10*plotB + 8);
plot(s,JAC), hold on,
xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
ylabel('JAC','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
end

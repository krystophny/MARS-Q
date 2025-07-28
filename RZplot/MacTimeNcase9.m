% plot time evolution of RE avalanche and VDE
% from the output of MARS-Q runs with NCASE=9 
% Note that this is a stand-alone code,
% without coupling to other Mac* procedures

kdev    = 54; 

MacReadParam

kprint = 0;
knormJ = 1; %normalize Jpara by Jpara(t=0)
tsteps = [100 200 300 400 500 600];

Tmax = 500000000.55; %maximal time [ms] to plot MARS-Q output
Tmin = 0.0e+04;
kscale = 1;
ks = 2;           %1: plot minor radius as s
                  %2: plot minor radius as psip=s.^2

plotnp = 3;       %select time slices for plotting
                  %0: plot all time slices
                  %1: plot end time step
                  %2: plot start and end time steps
                  %3: plot 21 time slices
                  %4: plot time slices with physical time step of dtphys in [ms]

stime  = 1.3;     %minor radius for ploting time trace; use rational surface(s) if stime>1 

knew = 1;         %1: basic

if ks==1, SLAB = '\psi_p^{1/2}'; end
if ks==2, SLAB = '\psi_p'; end

% get re-normalisation factor
fac_xn = 1.0;
fac_tm = 1.0;
fac_b1 = 1.0;
fac_jp = 1.0;
if kscale==1
   fac_xn = R0*1e+3;        %[mm]
   fac_tm = TAUA*1e+3;      %[ms]
   fac_b1 = B0*1e+4;        %[Gauss]
   fac_jp = R0*B0/(4e-7*pi)/1e+6;  %[MA]
end

ax_xn = 'displacement';
ax_tm = 't/\tau_A';ax_b1 = 'field';
ax_jp = 'parallel current';
ax_dt = '{\Delta}t/\tau_A';
if kscale==1
   ax_xn = 'displacement [mm]';
   ax_tm = 'time [ms]';
   ax_b1 = 'field [G]';
   ax_jp = 'parallel current [MA]';
   ax_dt = '{\Delta}t [ms]';
end
if knormJ==1, ax_jp='normalized parallel current'; end

Tmax = Tmax/fac_tm;

% get the radial grid
data = load([SDIR 'PROFEQ.OUT']);
s1   = data(:,1);
w0   = data(:,6);  w0e = w0(end);
q1   = data(:,2);
d0   = data(:,5);
Te   = data(:,11);
s    = (s1(1:end-1)+s1(2:end))/2;  
w0   = (w0(1:end-1)+w0(2:end))/2;
d0   = (d0(1:end-1)+d0(2:end))/2;

% read the time evolution data from the MARS-Q output
%filename='TIMEEVOL_delta_0.5';  kc=1;
%filename='TIMEEVOL_delta_fixIp_0.5';  kc=1;
%filename='TIMEEVOL_AS_kappa_fixIp_1.7';  kc=1;
%filename='TIMEEVOL_NS_kappa_1.7';  kc=1;
%filename='TIMEEVOL_AS_kappa_1.7';  kc=1;
%filename='TIMEEVOL_IpFac_1.2';  kc=1;
%filename='TIMEEVOL_SEED_2';  kc=1;
%filename='TIMEEVOL_14';  kc=1;  %fixIp
%filename='TIMEEVOL_20n';  kc=1;
%filename='TIMEEVOL_14d';  kc=3;
%filename='TIMEEVOL_07e';  kc=1;
filename='TIMEEVOL.OUT';  kc=1;
%filename='TIMEEVOL_SAVE';  kc=1;
%filename='../Time/TIMEEVOL_11_fixIp';  kc=1;

data = load([SDIR filename]);

delt = data(:,1); 
tt = cumsum(delt);
II = find(tt<Tmax & tt>Tmin);
data = data(II,:); t=tt(II);

jre  = data(:,3) + data(:,4)*i;
jpa  = data(:,5) + data(:,6)*i;
epa  = data(:,7)./data(:,2);  %=averaged E_para/E_c
solx = data(:,9);
soly = data(:,10);
jrec  = data(:,12) + data(:,13)*i;
jpac  = data(:,14) + data(:,15)*i;
epac  = data(:,16)./data(:,2);  %=averaged E_para/E_c

% specify line color
Nc = 5;
SS = '-';
cc = zeros(Nc,3);
for k=1:Nc
    cc(k,1) = (k-1)/(Nc-1);
    cc(k,3) = (Nc-k)/(Nc-1);
    cc(k,2) = (k-1)*(Nc-k)/(Nc-1)^2;
end
cc = [1 0 0; 0 0 1; 0 0 0; 0 1 0; 1 0 1];


if knormJ==1, fac_jp=1/real(jpa(1)); end
%====================================================================
% plot time steps
%====================================================================
hf = figure(1);
semilogy(tt*fac_tm,delt*fac_tm,'-','Color',cc(kc,:),'LineWidth',3), hold on,
xlabel(ax_tm,'FontSize',18,'FontWeight','Bold')
ylabel(ax_dt,'FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

%====================================================================
% plot time vs step
%====================================================================
hf = figure(2);
semilogy([1:length(tt)],tt*fac_tm,'-','Color',cc(kc,:),'LineWidth',3), hold on,
xlabel('time step','FontSize',18,'FontWeight','Bold')
ylabel(ax_tm,'FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

%====================================================================
% plot time traces of RE vs total parallel currents
%====================================================================
hf = figure(3);
semilogy(tt*fac_tm,real(jpa)*fac_jp,'-','Color',cc(kc,:),'LineWidth',3), hold on,
semilogy(tt*fac_tm,real(jre)*fac_jp,'--','Color',cc(kc,:),'LineWidth',3), hold on,
JRE_SEED_FRAC=real(jre(1))*fac_jp
%plot(tt*fac_tm,real(jpac)*fac_jp,'-.','Color',cc(kc,:),'LineWidth',3), hold on,
%plot(tt*fac_tm,real(jrec)*fac_jp,':','Color',cc(kc,:),'LineWidth',3), hold on,
xlabel(ax_tm,'FontSize',18,'FontWeight','Bold')
ylabel(ax_jp,'FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')


%====================================================================
% plot time traces of RE vs total parallel currents
% imaginery part just for checking
%====================================================================
if 1==0
hf = figure(4);
plot(tt*fac_tm,imag(jpa)*fac_jp,'-','Color',cc(kc,:),'LineWidth',3), hold on,
plot(tt*fac_tm,imag(jre)*fac_jp,'--','Color',cc(kc,:),'LineWidth',3), hold on,
xlabel(ax_tm,'FontSize',18,'FontWeight','Bold')
ylabel('Im(current)','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
end

%====================================================================
% plot time trace of averaged parallel electric field E_para/E_c
%====================================================================
hf = figure(5);
semilogy(tt*fac_tm,epa,'-','Color',cc(kc,:),'LineWidth',3), hold on,
%plot(tt*fac_tm,epac,'-.','Color',cc(kc,:),'LineWidth',3), hold on,
%plot([tt(1) tt(end)]*fac_tm,[0 0],'k-'), hold on
xlabel(ax_tm,'FontSize',18,'FontWeight','Bold')
ylabel('E_{||}/E_c','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')


%====================================================================
% plot time traces of solx and soly
%====================================================================
hf = figure(6);
semilogy(tt*fac_tm,solx,'-','Color',cc(kc,:),'LineWidth',3), hold on,
semilogy(tt*fac_tm,soly,'--','Color',cc(kc,:),'LineWidth',3), hold on,
xlabel(ax_tm,'FontSize',18,'FontWeight','Bold')
ylabel('max|X|,max|Y|','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
if length(tsteps)>0
   a=axis;
   for k=1:length(tsteps)
       semilogy([tt(tsteps(k)) tt(tsteps(k))]*fac_tm,[a(3) a(4)],'k--')
   end
end



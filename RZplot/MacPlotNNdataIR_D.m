% plot data from NN database for 3D equilibria

%SDIR = '/cscratch/liuy/WorkIR_D/Database/ITER/131025.0005/3D/n2/M/';
%SDIR = '/cscratch/liuy/WorkIR_D/Database/ITER/101007.0002/3D/n2/U/';
SDIR = '/cscratch/liuy/WorkIR_D/Database/DIII-D/157376/3D/n3/U/';
%SDIR = '/cscratch/liuy/WorkIR_D/Database/MAST/25075/3D/n4/U/';

kcheck = 0;

d  = load([SDIR 'NNdata_X1U.dat']);
NIN = round(d(1));
NOU = round(d(2));
NRZ = round(d(3));
NPR = round(d(4));
NBV = round(d(5));
NSS = round(d(6));
NUS = round(d(7));
NVM = round(d(8));

s1  = linspace(0,1,NPR-1);  s1=s1.^(1./(2+4*s1)); s1=[0 s1(2)/2 s1(2:end)];
s2  = linspace(0,1,NUS-1);  s2=s2.^(1./(2+4*s2)); s2=[0 s2(2)/2 s2(2:end)];
t9  = linspace(-pi,pi,NBV+1); t9=t9(1:end-1);

%plot plasma boundary
II = 9:NRZ+8;
RZ = d(II);
R1 = RZ(1:end-1);
Z1 = RZ(1:end-1);
t1 = linspace(-pi,pi,NRZ); t1=t1(:);
II = [1:8 10:24 26:NRZ-1];
Z1(II) = RZ(end) + (R1(II)-1).*tan(t1(II));
R1([9 25]) = 1;
hf=figure(1);
plot(R1,Z1,'bo','MarkerSize',9), hold on,
xlabel('R/R_0','FontSize',18,'FontWeight','Bold')
ylabel('Z/R_0','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
axis equal

%plot q-profile
II = NRZ+9:NRZ+NPR+7;
q  = d(II);
hf=figure(2);
plot(s1(1:end-1),q,'b--o','MarkerSize',9), hold on,
xlabel('s','FontSize',18,'FontWeight','Bold')
ylabel('q','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

%plot P-profile
II = NRZ+NPR+8:NRZ+2*NPR+7;
P  = d(II);
hf=figure(3);
plot(s1,P,'b--o','MarkerSize',9), hold on,
xlabel('s','FontSize',18,'FontWeight','Bold')
ylabel('P','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

%plot density-profile
II = NRZ+2*NPR+8:NRZ+3*NPR+7;
rho= d(II);
hf=figure(4);
plot(s1,rho,'b--o','MarkerSize',9), hold on,
xlabel('s','FontSize',18,'FontWeight','Bold')
ylabel('density','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

%plot rotation-profile
II = NRZ+3*NPR+8:NRZ+4*NPR+7;
w  = d(II);
hf=figure(5);
plot(s1,w,'b--o','MarkerSize',9), hold on,
xlabel('s','FontSize',18,'FontWeight','Bold')
ylabel('rotation','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

%plot resistivity-profile
II = NRZ+4*NPR+8:NRZ+5*NPR+6;
eta= d(II);
hf=figure(6);
plot(s1(1:end-1),eta,'b--o','MarkerSize',9), hold on,
xlabel('s','FontSize',18,'FontWeight','Bold')
ylabel('resistivity','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

%plot dBvac
II = NRZ+5*NPR+7:NRZ+5*NPR+NBV+6;
BvR= d(II);
II = NRZ+5*NPR+NBV+7:NRZ+5*NPR+2*NBV+6;
BvI= d(II);
hf=figure(7);
plot(t9,BvR,'r-',t9,BvI,'b--','MarkerSize',9), hold on,
xlabel('s','FontSize',18,'FontWeight','Bold')
ylabel('{\delta}B^n_{vac}','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
legend('Re','Im')

if kcheck==1
dd  = load([SDIR '../PROFEQ.OUT']);
ss = dd(:,1);
N  = length(ss);
dd  = load([SDIR 'SVD_X1U5.txt']);
S0  = dd(2,1:5);
UR0 = dd(3:N+2,1:5);
UI0 = dd(3:N+2,6:10);
VR0 = dd(N+3:end,1:5);
VI0 = dd(N+3:end,6:10);
end

%plot S
II = NIN+9:NIN+NSS+8;
S  = d(II);
hf=figure(8);
plot(S,'b--o','MarkerSize',9), hold on,
xlabel('SVD mode number','FontSize',18,'FontWeight','Bold')
ylabel('S-value','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

if kcheck==1
plot(S0,'b--+','MarkerSize',9), hold on,
end

%plot U
II = NIN+NSS+9:NIN+NSS+NSS*NUS*2+8;
U  = d(II);
UR = reshape(U(1:NSS*NUS),NUS,NSS);
UI = reshape(U(NSS*NUS+1:end),NUS,NSS);
hf=figure(9);
plot(s2,UR,'r-','LineWidth',2), hold on,
plot(s2,UI,'b--','LineWidth',2), hold on,
xlabel('s','FontSize',18,'FontWeight','Bold')
ylabel('U-vector','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

if kcheck==1
plot(ss,UR0,'g-','LineWidth',2), hold on,
plot(ss,UI0,'k--','LineWidth',2), hold on,
end

%plot V
II = NIN+NSS+NSS*NUS*2+9:NIN+NOU+8;
V  = d(II);
mm = -33:33;
VR = reshape(V(1:NSS*NVM),NVM,NSS);
VI = reshape(V(NSS*NVM+1:end),NVM,NSS);
hf=figure(10);
plot(mm,VR,'r-','LineWidth',2), hold on,
plot(mm,VI,'b--','LineWidth',2), hold on,
xlabel('m','FontSize',18,'FontWeight','Bold')
ylabel('V-vector','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

if kcheck==1
plot(mm,VR0,'g-','LineWidth',2), hold on,
plot(mm,VI0,'k--','LineWidth',2), hold on,
end

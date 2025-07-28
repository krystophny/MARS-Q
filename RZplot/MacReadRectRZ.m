% read data BPLASMA_RECTRZ and plot 

%SDIR = '/home/liuy/Work/RE-D3D/Data_184602/'; BNORM=1.6071e-4; R00=1.35247924E+00; Z00=1.45113933E-02;
SDIR = '/home/liuy/Work/RE-D3D/Data_185604/'; BNORM=1.3452e-4; R00=1.37963131E+00; Z00=-2.61982489E-02;

LCR     = 'b';
SLW     = 1; %LineWidth
KASCOT  = 0;
kplot2d = 1;
NZ     = 321;
NT      = 101;

plotB   = 0;   %1D plots
plotBB  = 3;   %2D plots
Rplot   = 1.35;   
Zplot   = 0.0;
n       = 1;
c0      = 0/6;  
WR      = 1.0;
WZ      = 0.0; 
WP      = 1-WR-WZ;
kRunA   = 0;  %0: no A-field involved 
              %1: read A from MARS-Q output file APLASMA_RECTRZ*, compute B from A
              %2: recompute A from B-field
kRunB   = 1;  %1: read B from MARS-Q output file BPLASMA_RECTRZ*
              %2: read B from MARS-Q output file BPLASMA_PSICHI* and BPLASMA_RECTRZ*
              %   and remapping B-field from PSI-CHI to R-Z mesh
              %3: plot BPLASMA in (s,theta) mesh
kSmAP   = 0;  %1: smooth A_PHI  

if kRunA==1 
   dataA = load([SDIR 'APLASMA_RECTRZ_C_n1_efitgrid']);
end;

if kRunB==1
   dataB = load([SDIR 'BPLASMA_RECTRZ.OUT']);
elseif kRunB==2
   dataB1 = load([SDIR 'BPLASMA_RECTRZ0001.OUT']);
   dataB2 = load([SDIR 'BPLASMA_PSICHI0001.OUT']);
elseif kRunB==3
   dataB = load([SDIR 'BPLASMA_PSICHI_PLS01']);
end

if kRunA==1
R    = dataA(:,1);
Z    = dataA(:,2);
AR   = dataA(:,3)+dataA(:,4)*i;
AZ   = dataA(:,5)+dataA(:,6)*i;
AP   = dataA(:,7)+dataA(:,8)*i;

II = find(R>0.1);
R    = R(II);
Z    = Z(II);
AR   = AR(II);
AZ   = AZ(II);
AP   = AP(II);

NR   = round(length(R)/NZ);

R    = reshape(R,NZ,NR); 
Z    = reshape(Z,NZ,NR);
AR   = reshape(AR,NZ,NR);
AZ   = reshape(AZ,NZ,NR);
AP   = reshape(AP,NZ,NR);
end

if kRunB==1 
R    = dataB(:,1);
Z    = dataB(:,2);
BR   = (dataB(:,3)+dataB(:,4)*i)*BNORM;
BZ   = (dataB(:,5)+dataB(:,6)*i)*BNORM;
BP   = (dataB(:,7)+dataB(:,8)*i)*BNORM;

II = find(R>0.1);
R    = R(II);
Z    = Z(II);
BR   = BR(II);
BZ   = BZ(II);
BP   = BP(II);

NR   = round(length(R)/NZ);

R    = reshape(R,NZ,NR); 
Z    = reshape(Z,NZ,NR);
BR   = reshape(BR,NZ,NR);
BZ   = reshape(BZ,NZ,NR);
BP   = reshape(BP,NZ,NR);
end

if kRunB==2
R2    = dataB2(:,1);
Z2    = dataB2(:,2);
BR2   = dataB2(:,3)+dataB2(:,4)*i;
BZ2   = dataB2(:,5)+dataB2(:,6)*i;
BP2   = dataB2(:,7)+dataB2(:,8)*i;

NS   = round(length(R2)/NT);

R2    = reshape(R2,NT,NS); 
Z2    = reshape(Z2,NT,NS);
BR2   = reshape(BR2,NT,NS);
BZ2   = reshape(BZ2,NT,NS);
BP2   = reshape(BP2,NT,NS);

R    = dataB1(:,1);
Z    = dataB1(:,2);
BR   = dataB1(:,3)+dataB1(:,4)*i;
BZ   = dataB1(:,5)+dataB1(:,6)*i;
BP   = dataB1(:,7)+dataB1(:,8)*i;

II = find(R>0.1);
R    = R(II);
Z    = Z(II);
BR   = BR(II);
BZ   = BZ(II);
BP   = BP(II);

NR   = round(length(R)/NZ);

R    = reshape(R,NZ,NR); 
Z    = reshape(Z,NZ,NR);
BR   = reshape(BR,NZ,NR);
BZ   = reshape(BZ,NZ,NR);
BP   = reshape(BP,NZ,NR);

BR = griddata(R2,Z2,BR2,R,Z);
BZ = griddata(R2,Z2,BZ2,R,Z);
BP = griddata(R2,Z2,BP2,R,Z);
end

if kRunB==3 
R    = dataB(:,1);   %=s grid
Z    = dataB(:,2);   %=theta grid
BR   = dataB(:,3)+dataB(:,4)*i;
BZ   = dataB(:,5)+dataB(:,6)*i;
BP   = dataB(:,7)+dataB(:,8)*i;

II = find(R>=0.94 & R<=1.06);
R = R(II);
Z = Z(II);
BR = BR(II);
BZ = BZ(II);
BP = BP(II);

NZ   = NT;
NR   = round(length(R)/NZ);

R    = reshape(R,NZ,NR); 
Z    = reshape(Z,NZ,NR);
BR   = reshape(BR,NZ,NR);
BZ   = reshape(BZ,NZ,NR);
BP   = reshape(BP,NZ,NR);

end

hR = 2*(R(1,2)-R(1,1));  hZ=2*(Z(2,1)-Z(1,1));  
IIR=[2:NR-1];            IIZ=[2:NZ-1];

if kRunA==2
% reconstruct A from B as weighted mixture of guage conditions AR=0, AZ=0, and AP=0
AR   = WR*BR*0 - WZ*cumsum(BP,1)*hZ/2 + WP*R.*BZ/(n*i);
AZ   = WR*cumsum(BP,2)*hR/2 + WZ*BZ*0 - WP*R.*BR/(n*i);
AP   =-WR*cumsum(R.*BZ,2)*hR/2./R + WZ*cumsum(BR,1)*hZ/2 + WP*BP*0;
end

%smooth A_PHI
if kSmAP>0 & kRunA>0
for k=1:20
   AP(IIZ,:) = c0*AP(IIZ+1,:)+c0*AP(IIZ-1,:)+(1-2*c0)*AP(IIZ,:);
end
end

% compute B from A
if kRunA>0
dAPdZ = AP;
dAPdZ(IIZ,IIR) =        c0*(AP(IIZ+1,IIR-1)-AP(IIZ-1,IIR-1))/hZ + ...
                                c0*(AP(IIZ+1,IIR+1)-AP(IIZ-1,IIR+1))/hZ + ...
                        (1-2*c0)*(AP(IIZ+1,IIR)-AP(IIZ-1,IIR))/hZ;
dAPdZ(1,:) = dAPdZ(2,:); dAPdZ(:,1) = dAPdZ(:,2);
dAPdZ(NZ,:) = dAPdZ(NZ-1,:); dAPdZ(:,NR) = dAPdZ(:,NR-1);
RAP = R.*AP; dRAPdR = RAP;
dRAPdR(IIZ,IIR) =         c0*(RAP(IIZ-1,IIR+1)-RAP(IIZ-1,IIR-1))/hR + ...
                                 c0*(RAP(IIZ+1,IIR+1)-RAP(IIZ+1,IIR-1))/hR + ...
                         (1-2*c0)*(RAP(IIZ,IIR+1)-RAP(IIZ,IIR-1))/hR;
dRAPdR(1,:) = dRAPdR(2,:); dRAPdR(:,1) = dRAPdR(:,2);
dRAPdR(NZ,:) = dRAPdR(NZ-1,:); dRAPdR(:,NR) = dRAPdR(:,NR-1);
dAZdR = AZ;
dAZdR(IIZ,IIR) =         c0*(AZ(IIZ-1,IIR+1)-AZ(IIZ-1,IIR-1))/hR + ...
                                 c0*(AZ(IIZ+1,IIR+1)-AZ(IIZ+1,IIR-1))/hR + ...
                         (1-2*c0)*(AZ(IIZ,IIR+1)-AZ(IIZ,IIR-1))/hR;
dAZdR(1,:) = dAZdR(2,:); dAZdR(:,1) = dAZdR(:,2);
dAZdR(NZ,:) = dAZdR(NZ-1,:); dAZdR(:,NR) = dAZdR(:,NR-1);
dARdZ = AP;
dARdZ(IIZ,IIR) =        c0*(AR(IIZ+1,IIR-1)-AR(IIZ-1,IIR-1))/hZ + ...
                                c0*(AR(IIZ+1,IIR+1)-AR(IIZ-1,IIR+1))/hZ + ...
                        (1-2*c0)*(AR(IIZ+1,IIR)-AR(IIZ-1,IIR))/hZ;
dARdZ(1,:) = dARdZ(2,:); dARdZ(:,1) = dARdZ(:,2);
dARdZ(NZ,:) = dARdZ(NZ-1,:); dARdZ(:,NR) = dARdZ(:,NR-1);

BR  =-dAPdZ + i*n*AZ./R;
BZ  = (dRAPdR-i*n*AR)./R;
BP  =-dAZdR + dARdZ;  
end

B    = sqrt(abs(BR).^2 + abs(BZ).^2 + abs(BP).^2);

DIVB = B;
RBR = R.*BR;
DIVB(IIZ,IIR) =  (c0*(RBR(IIZ-1,IIR+1)-RBR(IIZ-1,IIR-1))/hR + ...
                         c0*(RBR(IIZ+1,IIR+1)-RBR(IIZ+1,IIR-1))/hR + ...
                 (1-2*c0)*(RBR(IIZ,IIR+1)-RBR(IIZ,IIR-1))/hR)./R(IIZ,IIR) + ...
                         c0*(BZ(IIZ+1,IIR-1)-BZ(IIZ-1,IIR-1))/hZ + ...
                         c0*(BZ(IIZ+1,IIR+1)-BZ(IIZ-1,IIR+1))/hZ + ...
                 (1-2*c0)*(BZ(IIZ+1,IIR)-BZ(IIZ-1,IIR))/hZ+ ...
                         i*n*BP(IIZ,IIR)./R(IIZ,IIR);            
DIVB(2,:) = DIVB(3,:); DIVB(:,2) = DIVB(:,3);
DIVB(NZ-1,:) = DIVB(NZ-2,:); DIVB(:,NR-1) = DIVB(:,NR-2);
DIVB(1,:) = DIVB(3,:); DIVB(:,1) = DIVB(:,3);
DIVB(NZ,:) = DIVB(NZ-2,:); DIVB(:,NR) = DIVB(:,NR-2);

%compute a mask to set all fields outside surfs to be zero
Bmask = ones(size(R));
if 1==0
surfs = load([SDIR 'MacSurfS_D3D161015.1790']);
R00 = 1.6995;
surfs = surfs(1:end-1,:);
tt = atan2(Z,R-R00);
aa = sqrt((R-R00).^2 + Z.^2);
t  = atan2(surfs(:,2),surfs(:,1)-R00);
[t,I]=sort(t);  surfs = surfs(I,:);
a = sqrt((surfs(:,1)-R00).^2 + surfs(:,2).^2);

for k=1:NR
for j=1:NZ
    a1 = spline(t,a,tt(j,k));
    if (a1<aa(j,k)), Bmask(j,k)=0; end
end
end
end

if 1==1
   [I,J] = find(((R-R00).^2+(Z-Z00).^2)<0.00005);
   for k=1:length(I)
       Bmask(I(k),J(k)) = 0;
   end
end


% plot |B| in the R-Z plane
if kplot2d == 1
     hf=figure(10*plotBB + 0);
     pcolor(R,Z,B.*Bmask*1e+4), hold on, shading interp
     axis equal
     xlabel('R [m]','FontSize',14)
     ylabel('Z [m]','FontSize',14)
     title('|B| (Gauss)','FontSize',14)
     colorbar,  colormap(jet)
     ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
     maxAmpB = max(max(B.*Bmask))*1e+4
end

% save data for ASCOT
%note: need to change sign of BP to match ASCOT definition
if KASCOT==0 & 1==1
    %BP = -BP;
    Bn.r = R(1,:);
    Bn.z=  Z(:,1);
    Bn.Br = BR.*Bmask;
    Bn.Bz = BZ.*Bmask;
    Bn.Bphi =BP.*Bmask;
    save Bn Bn
end

%read data from ASCOT
if KASCOT==1
    eval(['load ' SSAV 'BnVAC.mat']);
    BR = Bn.Br;
    BZ = Bn.Bz;
    BP = Bn.Bphi;
end
if KASCOT==2
    eval(['load ' SSAV 'BnTOT.mat']);
    BR = Bn.Br;
    BZ = Bn.Bz;
    BP = Bn.Bphi;
end


if kplot2d == 1

% plot BR in the R-Z plane
hf=figure(10*plotBB + 1);
pcolor(R,Z,real(BR)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Re[\delta{B}_R]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

hf=figure(10*plotBB + 2);
pcolor(R,Z,imag(BR)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Im[\delta{B}_R]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

% plot BZ in the R-Z plane
hf=figure(10*plotBB + 3);
pcolor(R,Z,real(BZ)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Re[\delta{B}_Z]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

hf=figure(10*plotBB + 4);
pcolor(R,Z,imag(BZ)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Im[\delta{B}_Z]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

% plot Re(BP) in the R-Z plane
hf=figure(10*plotBB + 5);
pcolor(R,Z,real(BP.*Bmask)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Re[\delta{B}_\phi]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

hf=figure(10*plotBB + 6);
pcolor(R,Z,imag(BP.*Bmask)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Im[\delta{B}_\phi]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

end

BB = BR;

% plot BB at constant Z=Zplot
[Zmin,Imin] = min(abs(Z(:,1)-Zplot));
hf=figure(10*plotB + 4);
plot(R(Imin,:),real(BB(Imin,:)),[LCR '-'],R(Imin,:),imag(BB(Imin,:)),[LCR '--'],'LineWidth',SLW), hold on,
plot(R(Imin,:),abs(BB(Imin,:)),[LCR '-'],'LineWidth',3), hold on,
xlabel('R [m]','FontSize',16)
ylabel('\delta{B}_R','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend(['Real, Z=' num2str(Z(Imin,1)) 'm'],'Imag','abs')
a = axis;
%plot([Rplot Rplot],[a(3) a(4)],'k--'), hold on
plot([a(1) a(2)],[0 0],'k--'), hold on,

% plot BB at constant R=Rplot
[Rmin,Imin] = min(abs(R(1,:)-Rplot));
hf=figure(10*plotB + 5);
plot(Z(:,Imin),real(BB(:,Imin)),[LCR '-'],Z(:,Imin),imag(BB(:,Imin)),[LCR '--'],'LineWidth',SLW), hold on,
xlabel('Z [m]','FontSize',16)
ylabel('\delta{B}_R','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend(['Real, R=' num2str(R(1,Imin)) 'm'],'Imag')
a = axis;
%plot([Zplot Zplot],[a(3) a(4)],'k--'), hold on
plot([a(1) a(2)],[0 0],'k--'), hold on,

BB = BZ;

% plot BB at constant Z=Zplot
[Zmin,Imin] = min(abs(Z(:,1)-Zplot));
hf=figure(10*plotB + 6);
plot(R(Imin,:),real(BB(Imin,:)),[LCR '-'],R(Imin,:),imag(BB(Imin,:)),[LCR '--'],'LineWidth',SLW), hold on,
xlabel('R [m]','FontSize',16)
ylabel('\delta{B}_Z','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend(['Real, Z=' num2str(Z(Imin,1)) 'm'],'Imag')
a = axis;
%plot([Rplot Rplot],[a(3) a(4)],'k--'), hold on
plot([a(1) a(2)],[0 0],'k--'), hold on,

% plot BB at constant R=Rplot
[Rmin,Imin] = min(abs(R(1,:)-Rplot));
hf=figure(10*plotB + 7);
plot(Z(:,Imin),real(BB(:,Imin)),[LCR '-'],Z(:,Imin),imag(BB(:,Imin)),[LCR '--'],'LineWidth',SLW), hold on,
xlabel('Z [m]','FontSize',16)
ylabel('\delta{B}_Z','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend(['Real, R=' num2str(R(1,Imin)) 'm'],'Imag')
a = axis;
%plot([Zplot Zplot],[a(3) a(4)],'k--'), hold on
plot([a(1) a(2)],[0 0],'k--'), hold on,

BB = BP;

% plot BB at constant Z=Zplot
[Zmin,Imin] = min(abs(Z(:,1)-Zplot));
hf=figure(10*plotB + 8);
plot(R(Imin,:),real(BB(Imin,:)).*Bmask(Imin,:),[LCR '-'],R(Imin,:),imag(BB(Imin,:)).*Bmask(Imin,:),[LCR '--'],'LineWidth',SLW), hold on,
xlabel('R [m]','FontSize',16)
ylabel('\delta{B}_\phi','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend(['Real, Z=' num2str(Z(Imin,1)) 'm'],'Imag')
a = axis;
%plot([Rplot Rplot],[a(3) a(4)],'k--'), hold on
plot([a(1) a(2)],[0 0],'k--'), hold on,

% plot BB at constant R=Rplot
[Rmin,Imin] = min(abs(R(1,:)-Rplot));
hf=figure(10*plotB + 9);
plot(Z(:,Imin),real(BB(:,Imin)).*Bmask(:,Imin),[LCR '-'],Z(:,Imin),imag(BB(:,Imin)).*Bmask(:,Imin),[LCR '--'],'LineWidth',SLW), hold on,
xlabel('Z [m]','FontSize',16)
ylabel('\delta{B}_\phi','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend(['Real, R=' num2str(R(1,Imin)) 'm'],'Imag')
a = axis;
%plot([Zplot Zplot],[a(3) a(4)],'k--'), hold on
plot([a(1) a(2)],[0 0],'k--'), hold on,

if kplot2d == 1
    
BB = DIVB;

% plot DIVB at constant Z=Zplot
[Zmin,Imin] = min(abs(Z(:,1)-Zplot));
hf=figure(10*plotB + 10);
plot(R(Imin,:),real(BB(Imin,:)),[LCR '-'],R(Imin,:),imag(BB(Imin,:)),[LCR '--'],'LineWidth',SLW), hold on,
xlabel('R [m]','FontSize',16)
ylabel('div(\delta{B})','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend(['Real, Z=' num2str(Z(Imin,1)) 'm'],'Imag')
a = axis;
%plot([Rplot Rplot],[a(3) a(4)],'k--'), hold on
plot([a(1) a(2)],[0 0],'k--'), hold on,

% plot DIVB at constant R=Rplot
[Rmin,Imin] = min(abs(R(1,:)-Rplot));
hf=figure(10*plotB + 11);
plot(Z(:,Imin),real(BB(:,Imin)),[LCR '-'],Z(:,Imin),imag(BB(:,Imin)),[LCR '--'],'LineWidth',SLW), hold on,
xlabel('Z [m]','FontSize',16)
ylabel('div(\delta{B})','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend(['Real, R=' num2str(R(1,Imin)) 'm'],'Imag')
a = axis;
%plot([Zplot Zplot],[a(3) a(4)],'k--'), hold on
plot([a(1) a(2)],[0 0],'k--'), hold on,

end

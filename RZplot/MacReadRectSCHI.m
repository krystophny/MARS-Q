% read data SCHIMESH_RECTRZ and plot 

SDIR   = '/home/liuy/Work/D3D161015/Data/';

LCR     = 'b';
SLW     = 1; %LineWidth
kplot1d = 0; 
kplot2d = 2;
NZ      = 341; 
Rplot   = 1.6995;   
Zplot   = 0.0;

data = load([SDIR 'SCHIMESH_RECTRZ_n1_new']);

Raxis = 1.74901357E+00;
Zaxis =-1.20923981E-01;

R    = data(:,1);
Z    = data(:,2);
S    = data(:,3);
C    = data(:,4);

NR   = round(length(R)/NZ);

R    = reshape(R,NZ,NR); 
Z    = reshape(Z,NZ,NR);
S    = reshape(S,NZ,NR);
C    = reshape(C,NZ,NR);

if kplot2d > 0

% plot S-mesh in the R-Z plane
hf=figure(10*kplot2d + 1);
pcolor(R,Z,S), hold on, shading interp, 
colorbar,  colormap(hot)
plot(Raxis,Zaxis,'w+','MarkerSize',12)
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('s-mesh','FontSize',14)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

% plot CHI-mesh in the R-Z plane
hf=figure(10*kplot2d + 2);
pcolor(R,Z,C), hold on, shading interp, 
colorbar,  colormap(hot)
plot(Raxis,Zaxis,'w+','MarkerSize',12)
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('{\chi}-mesh','FontSize',14)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

end

if kplot1d > 0

BB = S;

% plot BB at constant Z=Zplot
[Zmin,Imin] = min(abs(Z(:,1)-Zplot));
hf=figure(10*kplot1d + 1);
plot(R(Imin,:),BB(Imin,:),[LCR '-'],'LineWidth',SLW), hold on,
xlabel('R [m]','FontSize',16)
ylabel('s-mesh','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

% plot BB at constant R=Rplot
[Rmin,Imin] = min(abs(R(1,:)-Rplot));
hf=figure(10*kplot1d + 2);
plot(Z(:,Imin),BB(:,Imin),[LCR '-'],'LineWidth',SLW), hold on,
xlabel('Z [m]','FontSize',16)
ylabel('s-mesh','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

BB = C;

% plot BB at constant Z=Zplot
[Zmin,Imin] = min(abs(Z(:,1)-Zplot));
hf=figure(10*kplot1d + 3);
plot(R(Imin,:),BB(Imin,:),[LCR '-'],'LineWidth',SLW), hold on,
xlabel('R [m]','FontSize',16)
ylabel('{\chi}-mesh','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

% plot BB at constant R=Rplot
[Rmin,Imin] = min(abs(R(1,:)-Rplot));
hf=figure(10*kplot1d + 4);
plot(Z(:,Imin),BB(:,Imin),[LCR '-'],'LineWidth',SLW), hold on,
xlabel('Z [m]','FontSize',16)
ylabel('{\chi}-mesh','FontSize',16)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

end

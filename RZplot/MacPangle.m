% plot poloidal angles in EQAC(GEOM) and PEST coordinates

SDIR = '/.automount/funsrv1/root/home/yliu/MAST_ADD/';
m    = 2; %which surface to plot

% EQAC angle
data = load([SDIR 'PANGLE_EQAC.IN']);

M = data(1);  %number of rational surfaces
N = data(2);  %number of points at each surface

chi2 = data((m-1)*N+3:m*N+2);
chi1 = linspace(0,2*pi,N);

hf=figure(1);
plot(chi1,chi2,'b--','LineWidth',2), hold on,

% PEST angle
data = load([SDIR 'PANGLE_PEST.IN']);

M = data(1);  %number of rational surfaces
N = data(2);  %number of points at each surface

chi2 = data((m-1)*N+3:m*N+2);
chi1 = linspace(0,2*pi,N);

hf=figure(1);
plot(chi1,chi2,'r-','LineWidth',2), hold on,

xlabel('\chi_{MARS}','FontSize',18,'FontWeight','Bold'),
ylabel('\chi','FontSize',18,'FontWeight','Bold'),

legend('GEOM','PEST')

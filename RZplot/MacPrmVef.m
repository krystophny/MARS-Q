% compute the n=1 vacuum error field
% due to horizontal displacement of PF coil

function [Bnm,Btm]=MacPrmVef(RR,ZZ,IB,T,IIT)

global Mac

% specify source coil location, horizontal displacement
% coil PF current 
% for F7A coil in D3D141069
Rc = 2.375; %[m]
Zc = 1.124; %[m]
Dc = 0.01;  %[m] horizontal displacement
Ic = 56*2.55e+3; %[A-turn]

if 1==0
% plot top view of coil geometry
pp = linspace(-pi,pi,201);
cosp = cos(pp);
sinp = sin(pp);
hf=figure(20);
plot(1.070*cosp,1.070*sinp,'b-','LineWidth',1), hold on,
plot(2.276*cosp,2.276*sinp,'b-','LineWidth',1), hold on,
plot(Rc*cosp,Rc*sinp,'r-','LineWidth',3), hold on,
plot(Rc*cosp+Dc,Rc*sinp,'r--','LineWidth',3), hold on,
xlabel('X [m]','FontSize',16,'FontWeight','Bold')
ylabel('Y [m]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
axis equal
end

% get control surface where the normal component of the n=1 error field is computed
R = RR(IB,:)*Mac.R0EXP;
Z = ZZ(IB,:)*Mac.R0EXP;
R      = R(IIT);
Z      = Z(IIT);

% compute the n=1 error field components
% first compute HatQ_{0,1,2}
R22 = R.^2 + Rc^2 + (Z-Zc).^2;
kappa = 2*R*Rc./R22;
kappa_min = min(kappa);
kappa_max = max(kappa);
ka = linspace(kappa_min,kappa_max,101)';
a = linspace(0,pi,201); 
ha = a(2);
a = (a(1:end-1)+a(2:end))/2;
HQ0 = (1-ka*cos(a)).^(-5/2);
HQ1 = ones(size(ka))*cos(a).*HQ0;
HQ2 = ones(size(ka))*cos(2*a).*HQ0;
HQ0 = sum(HQ0,2)*ha;
HQ1 = sum(HQ1,2)*ha;
HQ2 = sum(HQ2,2)*ha;

HQ0 = spline(ka,HQ0,kappa);
HQ1 = spline(ka,HQ1,kappa);
HQ2 = spline(ka,HQ2,kappa);

% compute BRn1,BZn1 in Gauss
mu0 = 4e-3*pi;
R25 = R22.^(5/2);

BRn1 = 3*mu0*Ic*Dc/8/pi*Rc*(Z-Zc)./R25.*(2*R.*HQ1-Rc*HQ0-Rc*HQ2);
BZn1 = -mu0*Ic*Dc/8/pi*Rc./R25.*(2*(3*R.^2+3*Rc^2-R22).*HQ1 - 7*Rc*R.*HQ0 - Rc*R.*HQ2);

% compute normal component Bn and tangential component Bt
% first compute dR/dT, and dZ/dT
dRT = diff(R)./diff(T);
x = (T(2:end)+T(1:end-1))/2;
x = [x(end)-2*pi x x(1)+2*pi];
dRT = [dRT(end) dRT dRT(1)];
dR = pchip(x,dRT,T);
dZT = diff(Z)./diff(T);
dZT = [dZT(end) dZT dZT(1)];
dZ = pchip(x,dZT,T);

%compute physical Bn & Bt from BR & BZ
nA = sqrt(dR.^2+dZ.^2);
Bn = (BRn1.*dZ - BZn1.*dR)./nA;
Bt = (BRn1.*dR + BZn1.*dZ)./nA;

%Fourier decompose Bn and Bt along T angle
%get Gauss quadrature points for Fourier decomposition
[z,w] = MacGaussQuad1D(2);
x0 = (T(1:end-1)+T(2:end))*0.5;
h2 = diff(T)*0.5;
xx = z'*h2 + ones(size(z'))*x0;
wh = w'*h2;
xx = xx(:);  wh=wh(:);

expmt = exp(-i*xx*Mac.Mm')/(2*pi);
Bnx   = spline(T,Bn,xx').*wh';
Btx   = spline(T,Bt,xx').*wh';
Bnm = Bnx*expmt; 
Btm = Btx*expmt; 
Bnm = Bnm(:);
Btm = Btm(:);





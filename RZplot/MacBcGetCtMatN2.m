% a new and better version for computing coordinate transform matrices
%
function [A1] = MacBcGetCtMatN2(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian)

global Mac Fc

%find minor radius for coupling surface
[smin,II] = min(abs(Mac.s-Fc.rs));
coupling_surface_radius = Mac.s(II)

%get quantities on coupling surface
Cw = Mac.chi;
Rw = R(II,:)*Mac.R0EXP;
Zw = Z(II,:)*Mac.R0EXP;

%find geometrical origin to define theta angle
Zmin=min(Zw);  Zmax=max(Zw); Z0=0.5*(Zmin+Zmax);
Rmin=min(Rw);  Rmax=max(Rw); R0=0.5*(Rmin+Rmax);

%save coupling surfaces
Rw2 = R(II-1,:)*Mac.R0EXP;
Zw2 = Z(II-1,:)*Mac.R0EXP;
Rw2 = (Rw+Rw2)/2;
Zw2 = (Zw+Zw2)/2;
res = [Rw' Zw' Rw2' Zw2'];
save MacDataBcS res -ascii -double
eval(['!mv MacDataBcS ' Fc.DIR])

%compute geometrical angle
Tw = atan2(Zw-Z0,Rw-R0);
I = find(diff(Tw)<-pi); 
if length(I)==1
  Tw(I+1:end) = Tw(I+1:end) + 2*pi;  
end

%compute dTheta/dChi
x0 = (Cw(1:end-1)+Cw(2:end))/2;
y0 = diff(Tw)./diff(Cw);
dT = spline(x0,y0,Cw);

%compute C-factor by computing first (dRw/dTw,dZw/dTw)
%get Gauss quadrature points for Fourier decomposition
[z,w] = MacGaussQuad1D(4);
x0 = (Tw(1:end-1)+Tw(2:end))*0.5;
h2 = diff(Tw)*0.5;
xx = z'*h2 + ones(size(z'))*x0;
wh = w'*h2;
xx = xx(:);  wh=wh(:);

%Fourier decompose (Rw,Zw) in Tw-angle
expmt = exp(-i*xx*[-Mac.Nm2:Mac.Nm2])/(2*pi);
Rx    = spline(Tw,Rw,xx').*wh';
Zx    = spline(Tw,Zw,xx').*wh';
Rm = Rx*expmt; 
Zm = Zx*expmt;

%compute (dRw/dTw,dZw/dTw) via Fourier harmonics
expmt = exp(i*[-Mac.Nm2:Mac.Nm2]'*Tw);
dR = i*([-Mac.Nm2:Mac.Nm2].*Rm)*expmt;
dZ = i*([-Mac.Nm2:Mac.Nm2].*Zm)*expmt;
dR = real(dR);
dZ = real(dZ);

%compute C-factor
C = sqrt(dR.^2+dZ.^2);

%compute coordinate transform matrices
%using Gauss quadrature rules for integration
x0 = (Cw(1:end-1)+Cw(2:end))*0.5;
h2 = diff(Cw)*0.5;
xx = z'*h2 + ones(size(z'))*x0;
wh = w'*h2;
xx = xx(:);  wh=wh(:); wh=wh';
Tx    = spline(Cw,Tw,xx);
Rx    = spline(Cw,Rw,xx');
Cx    = spline(Cw,C,xx');
dTx   = spline(Cw,dT,xx');

mm = Fc.mm;  M=length(mm);
A1 = zeros(M,M);  
oo = ones(size(mm));

for m=1:M
  expmt = exp(i*Tx*mm-i*xx*(mm(m)*oo))/(2*pi);
  A1(m,:) = (wh.*Cx.*Rx.*dTx)*expmt;
end

%plotting ...
if Fc.plotC>0
  figure(1+10*Fc.plotC)
  plot(Cw,Tw),
  xlabel('\chi'), ylabel('\theta')

  figure(2+10*Fc.plotC)
  plot(Cw,Rw),
  xlabel('\chi'), ylabel('R')

  figure(3+10*Fc.plotC)
  plot(Cw,Zw),
  xlabel('\chi'), ylabel('Z')

  figure(5+10*Fc.plotC)
  plot(Cw,C),
  xlabel('\chi'), ylabel('C')
end


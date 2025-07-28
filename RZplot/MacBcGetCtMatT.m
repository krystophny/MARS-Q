function [A2a,A2b,A3] = MacBcGetCtMatT(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian)

global Mac Fc

%find minor radius for coupling surface
[smin,II] = min(abs(Mac.s-Fc.rs));

%get quantities on coupling surface
Cw = Mac.chi;
Rw = (R(II,:)+R(II-1,:))/2*Mac.R0EXP;
Zw = (Z(II,:)+Z(II-1,:))/2*Mac.R0EXP;
Jw = (jacobian(II,:)+jacobian(II-1,:))/2*Mac.R0EXP.^3;
G12w1 = (dRds(II,:).*dRdchi(II,:) + dZds(II,:).*dZdchi(II,:))*Mac.R0EXP.^2;
G22w1 = (dRdchi(II,:).^2 + dZdchi(II,:).^2)*Mac.R0EXP.^2;
G12w2 = (dRds(II-1,:).*dRdchi(II-1,:) + dZds(II-1,:).*dZdchi(II-1,:))*Mac.R0EXP.^2;
G22w2 = (dRdchi(II-1,:).^2 + dZdchi(II-1,:).^2)*Mac.R0EXP.^2;
G12w = (G12w1+G12w2)/2;
G22w = (G22w1+G22w2)/2;

%find geometrical origin to define theta angle
Zmin=min(Zw);  Zmax=max(Zw); Z0=0.5*(Zmin+Zmax)
Rmin=min(Rw);  Rmax=max(Rw); R0=0.5*(Rmin+Rmax)

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
Jx    = spline(Cw,Jw,xx');
G12x    = spline(Cw,G12w,xx');
G22x    = spline(Cw,G22w,xx');
Cx    = spline(Cw,C,xx');
dTx   = spline(Cw,dT,xx');

mm = Fc.mm;  M=length(mm);
A1 = zeros(M,M);  A2a=A1; A2b=A1; A3=A1;
oo = ones(size(mm));

for m=1:M
  expmt = exp(i*xx*mm-i*Tx*(mm(m)*oo))/(2*pi);
  A2a(m,:) = (wh.*G12x./Cx./Jx)*expmt;
  A2b(m,:) = (wh.*G22x./Cx./Jx)*expmt;
  A3(m,:) = (wh.*Rx./Jx.*dTx)*expmt;
end

%invert matrices A1,A2,A3
A2b = inv(A2b);
A3 = -inv(A3);

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

  figure(4+10*Fc.plotC)
  plot(Cw,Jw),
  xlabel('\chi'), ylabel('Jacobian')

  figure(5+10*Fc.plotC)
  plot(Cw,C),
  xlabel('\chi'), ylabel('C')
end


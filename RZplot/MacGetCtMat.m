function [A1,A2,A3,II] = MacBcGetCtMat(R,Z,dRdchi,dZdchi)

global Mac 

%find minor radius for coupling surface
[smin,II] = min(abs(Mac.s-Mac.rbn));

%get quantities on coupling surface
Cw = Mac.chi;
Rw = R(II,:);
Zw = Z(II,:);
G22w = sqrt(dRdchi(II,:).^2 + dZdchi(II,:).^2);

%find geometrical origin to define theta angle
%Zmin=min(Zw);  Zmax=max(Zw); Z0=0.5*(Zmin+Zmax)
%Rmin=min(Rw);  Rmax=max(Rw); R0=0.5*(Rmin+Rmax)
R0 = R(1,1);
Z0 = Z(1,1);

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

kdRZ = 2;
if kdRZ == 1
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
else
%use simple finite difference to compute dR,dZ
T0 = (Tw(1:end-1)+Tw(2:end))/2;
dR = diff(Rw)./diff(Tw);
dZ = diff(Zw)./diff(Tw);
dR = spline(T0,dR,Tw);
dZ = spline(T0,dZ,Tw);
end

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
G22x    = spline(Cw,G22w,xx');
Cx    = spline(Cw,C,xx');
dTx   = spline(Cw,dT,xx');

mm = Mac.Mm';  M=length(mm);
A1 = zeros(M,M);  A2=A1; A3=A1;
oo = ones(size(mm));

for m=1:M
  expmt = exp(i*xx*mm-i*Tx*(mm(m)*oo))/(2*pi);
  A1(m,:) = (wh./Cx./Rx)*expmt;
  expmt = exp(i*xx*(mm-mm(m)))/(2*pi);
  A2(m,:) = (wh./G22x./Rx)*expmt;
  expmt = exp(i*Tx*mm-i*xx*(mm(m)*oo))/(2*pi);
  A3(m,:) = (wh.*Cx.*Rx.*dTx)*expmt;
end

if 1==0
  figure(30+0)
  plot(Cw,Tw,Mac.SS), hold on,
  xlabel('\chi'), ylabel('\theta')


  figure(30+1)
  plot(Tw,Rw,Mac.SS), hold on,
  xlabel('\theta'), ylabel('R')

  figure(30+2)
  plot(Tw,Zw,Mac.SS), hold on,
  xlabel('\theta'), ylabel('Z')

  figure(30+3)
  plot(Tw,C,Mac.SS), hold on,
  xlabel('\theta'), ylabel('C')
end

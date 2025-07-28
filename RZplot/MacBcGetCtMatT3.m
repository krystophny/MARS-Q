% a new and better version for computing coordinate transform matrices
%
function [A2a,A2b,A3] = MacBcGetCtMatT3(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian)

global Mac Fc

%find minor radius for coupling surface
[smin,II] = min(abs(Mac.s-Fc.rs));

%get quantities on coupling surface
Cw = Mac.chi;
Rw1 = R(II,:)*Mac.R0EXP;
Zw1 = Z(II,:)*Mac.R0EXP;
Rw2 = R(II-1,:)*Mac.R0EXP;
Zw2 = Z(II-1,:)*Mac.R0EXP;
Rw = (Rw1+Rw2)/2;
Zw = (Zw1+Zw2)/2;
Jw = (jacobian(II,:)+jacobian(II-1,:))/2*Mac.R0EXP.^3;
G12w1 = (dRds(II,:).*dRdchi(II,:) + dZds(II,:).*dZdchi(II,:))*Mac.R0EXP.^2;
G22w1 = (dRdchi(II,:).^2 + dZdchi(II,:).^2)*Mac.R0EXP.^2;
G12w2 = (dRds(II-1,:).*dRdchi(II-1,:) + dZds(II-1,:).*dZdchi(II-1,:))*Mac.R0EXP.^2;
G22w2 = (dRdchi(II-1,:).^2 + dZdchi(II-1,:).^2)*Mac.R0EXP.^2;
G12w = (G12w1+G12w2)/2;
G22w = (G22w1+G22w2)/2; G22w = sqrt(G22w);

%find geometrical origin to define theta angle
Zmin=min(Zw1);  Zmax=max(Zw1); Z0=0.5*(Zmin+Zmax);
Rmin=min(Rw1);  Rmax=max(Rw1); R0=0.5*(Rmin+Rmax);

%compute geometrical angle
Tw = atan2(Zw-Z0,Rw-R0);
I = find(diff(Tw)<-pi); 
if length(I)==1
  Tw(I+1:end) = Tw(I+1:end) + 2*pi;  
end

%compute coordinate transform matrices
%using Gauss quadrature rules for integration
[z,w] = MacGaussQuad1D(4);
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

mm = Fc.mm;  M=length(mm);
A1 = zeros(M,M);  A2a=A1; A2b=A1; A3=A1;
oo = ones(size(mm));

for m=1:M
  expmt = exp(i*Tx*mm-i*xx*(mm(m)*oo))/(2*pi);
  A2a(m,:) = (wh./G22x.*Jx)*expmt;
  A2b(m,:) = (wh.*G12x./G22x.*Rx)*expmt;
  A3(m,:) = -(wh./Rx.*Jx)*expmt;
end

%plotting ...
if Fc.plotC>0
  figure(1+10*Fc.plotC)
  plot(Cw,Tw,'b--'), hold on,
  xlabel('\chi'), ylabel('\theta')

  figure(2+10*Fc.plotC)
  plot(Cw,Rw,'b--'), hold on,
  xlabel('\chi'), ylabel('R')

  figure(3+10*Fc.plotC)
  plot(Cw,Zw,'b--'), hold on
  xlabel('\chi'), ylabel('Z')

  figure(5+10*Fc.plotC)
  plot(Cw,Jw,'b--'), hold on,
  xlabel('\chi'), ylabel('Jacobian')
end


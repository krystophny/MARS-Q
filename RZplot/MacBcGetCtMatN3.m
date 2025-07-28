% a new and better version for computing coordinate transform matrices
%
function [A1] = MacBcGetCtMatN3(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian)

global Mac Fc

CheckA1 = 0;

%find minor radius for coupling surface
[smin,II] = min(abs(Mac.s-Fc.rs));
coupling_surface_radius = Mac.s(II)

%get quantities on coupling surface
Cw = Mac.chi;
Rw = R(II,:)*Mac.R0EXP;
Zw = Z(II,:)*Mac.R0EXP;
G22w = (dRdchi(II,:).^2 + dZdchi(II,:).^2)*Mac.R0EXP.^2;
G22w = sqrt(G22w);

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
G22x  = spline(Cw,G22w,xx');
dTx   = spline(Cw,dT,xx');

mm = Fc.mm;  M=length(mm);
A1 = zeros(M,M);  
oo = ones(size(mm));

for m=1:M
  expmt = exp(i*xx*mm-i*Tx*(mm(m)*oo))/(2*pi);
  A1(m,:) = (wh./Rx./G22x.*dTx)*expmt;
end

if CheckA1 > 0
   A1inv = A1;  
   for m=1:M
       expmt = exp(i*Tx*mm-i*xx*(mm(m)*oo))/(2*pi);
       A1inv(m,:) = (wh.*Rx.*G22x)*expmt;
   end

   R = A1*A1inv;
   V = 1./real(diag(R));
   D = diag(V.^2);
   A1 = A1*D;

   R = A1*A1inv;
   figure(1)
   plot(mm,real(diag(R)),'r-o'), hold on,
   plot(mm,imag(diag(R)),'b-x'), hold on,
end


%plotting ...
if Fc.plotC>0
  figure(1+10*Fc.plotC)
  plot(Cw,Tw,'r-'), hold on,
  xlabel('\chi'), ylabel('\theta')

  figure(2+10*Fc.plotC)
  plot(Cw,Rw,'r-'), , hold on,
  xlabel('\chi'), ylabel('R')

  figure(3+10*Fc.plotC)
  plot(Cw,Zw,'r-'), , hold on,
  xlabel('\chi'), ylabel('Z')

  figure(4+10*Fc.plotC)
  plot(Cw,dT,'r-'), , hold on,
  xlabel('\chi'), ylabel('dT')
end


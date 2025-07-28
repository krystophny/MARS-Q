function [y] = MacGetDenPumpOut(R,Vn,jac,den)

global Mac

Iedge = Mac.Ns1-3;

VV = real(exp(i*Mac.n*Mac.phi)*Vn(Iedge,:));
II = find(VV<0); VV(II)=0;

Vchi = sum(VV,1)*(Mac.phi(2)-Mac.phi(1));
den1 = sum(Vchi.*jac(Iedge,:))*(Mac.chi(2)-Mac.chi(1))*den(Iedge);

h    = diff(Mac.s(1:Mac.Ns1));
jac1 = sum(jac,2)*(Mac.chi(2)-Mac.chi(1));
jac2 = (jac1(1:Mac.Ns1-1)+jac1(2:Mac.Ns1))/2;
d    = (den(1:Mac.Ns1-1)+den(2:Mac.Ns1))/2;
den2 = sum(d.*jac2.*h)*2*pi;

%get aspect ratio
Rmin = min(R(Mac.Ns1,:));
Rmax = max(R(Mac.Ns1,:));
A    = (Rmin+Rmax)/(Rmax-Rmin);

den2 = den2*Mac.R0EXP/A;

y    = den1/den2;

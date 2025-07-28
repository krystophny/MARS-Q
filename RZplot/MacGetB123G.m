function [CHIG,WG,B1G,B2G,B3G] = MacGetB123G(TC,BrC,BzC,BphiC,R,Z,dRds,dZds,dRdchi,dZdchi,jacobian)

%% compute B1,B2,B3 from BrC,BzC,BphiC, at Gaussian points 

global Mac

% Compute Gaussian quadrature
n = 4;  % number of quadrature points
[x,w] = MacGaussQuad1D(n);

chi = Mac.chi;
N = length(chi)-1;  h = chi(2)-chi(1);
CHIG = x(:)*ones(1,N)*h/2 + ones(n,1)*chi(1:N)+h/2;
WG   = w(:);
 
% spline interporations
[smin,II] = min(abs(Mac.s-Mac.rw(1)));
if smin<(Mac.s(II)-Mac.s(II-1))*0.1
  RO=R(II,:);
  ZO=Z(II,:);
  dRdsO=dRds(II,:);
  dZdsO=dZds(II,:);
  dRdchiO=dRdchi(II,:);
  dZdchiO=dZdchi(II,:);
  jacobianO=jacobian(II,:);
else
  RO=(R(II,:)+R(II-1,:))/2;
  ZO=(Z(II,:)+Z(II-1,:))/2;
  dRdsO=(dRds(II,:)+dRds(II-1,:))/2;
  dZdsO=(dZds(II,:)+dZds(II-1,:))/2;
  dRdchiO=(dRdchi(II,:)+dRdchi(II-1,:))/2;
  dZdchiO=(dZdchi(II,:)+dZdchi(II-1,:))/2;
  jacobianO=(jacobian(II,:)+jacobian(II-1,:))/2;
end

RG = spline(chi,RO,CHIG(:));
ZG = spline(chi,ZO,CHIG(:));
dRdsG = spline(chi,dRdsO,CHIG(:));
dZdsG = spline(chi,dZdsO,CHIG(:));
dRdchiG = spline(chi,dRdchiO,CHIG(:));
dZdchiG = spline(chi,dZdchiO,CHIG(:));
jacobianG = spline(chi,jacobianO,CHIG(:));
TG = atan2(ZG,RG-1.0);
BrG = spline(TC,BrC,TG);
BzG = spline(TC,BzC,TG);
BphiG = spline(TC,BphiC,TG);

% convert from (r,z,phi) to (1,2,3)
B1G = RG.*(dZdchiG.*BrG - dRdchiG.*BzG)/Mac.B0EXP;
B2G = RG.*( -dZdsG.*BrG + dRdsG.*BzG)/Mac.B0EXP;
B3G = jacobianG.*BphiG/./RG/Mac.B0EXP;

B1G = reshape(B1G,n,N);
B2G = reshape(B2G,n,N);
B3G = reshape(B3G,n,N);
 

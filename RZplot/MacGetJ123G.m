function [CHIG,WG,J1G,J2G,J3G] = MacGetJ123G(TC,JrC,JzC,JphiC,R,Z,dRds,dZds,dRdchi,dZdchi,jacobian)

%% compute J1,J2,J3 from JrC,JzC,JphiC, at Gaussian points 

global Mac

% Compute Gaussian quadrature
n = 4;  % number of quadrature points
[x,w] = MacGaussQuad1D(n);

chi = Mac.chi;
N = length(chi)-1;  h = chi(2)-chi(1);
CHIG = x(:)*ones(1,N)*h/2 + ones(n,1)*chi(1:N)+h/2;
WG   = w(:);
 
% spline interporations
[smin,II] = min(abs(Mac.s-Mac.rs(1)));
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
JrG = spline(TC,JrC,TG);
JzG = spline(TC,JzC,TG);
JphiG = spline(TC,JphiC,TG);

% convert from (r,z,phi) to (1,2,3)
J1G = RG.*(dZdchiG.*JrG - dRdchiG.*JzG);
J2G = RG.*( -dZdsG.*JrG + dRdsG.*JzG);
J3G = jacobianG.*JphiG./RG;

J1G = reshape(J1G,n,N);
J2G = reshape(J2G,n,N);
J3G = reshape(J3G,n,N);
 

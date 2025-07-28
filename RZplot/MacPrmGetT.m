% compute geometrical angle T at control surface

function [T,IIT]=MacPrmGetT(RR,ZZ,IB)

global Mac

% get control surface where the normal component of the n=1 error field is computed
R = RR(IB,:);
Z = ZZ(IB,:);

% compute geometrical poloidal angle T
Zmin=min(Z);  Zmax=max(Z); Z0=0.5*(Zmin+Zmax);
Rmin=min(R);  Rmax=max(R); R0=0.5*(Rmin+Rmax);

T      = atan2(Z-Z0,R-R0);
[T,IIT] = sort(T);


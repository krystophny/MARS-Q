function FEEDI=MacACADcoilcurrent(n,R0,B0)

global Acad

% coil current fC=|In/I0| for different RMP coil configurations 
% for DIII-D I-coils

% Generic expression for computing fA=I0/In 
% for a given array of coil currents
% assuming that all coils are equally distributed, and
% all coils have the same toroidal extension

if length(Acad)>0
if strcmp(Acad.DEVICE,'D3D')
dphi = 60/180*pi;  %D3D I-coil and C-coil
end
if strcmp(Acad.DEVICE,'MAST')
dphi = 0.4061;     %MAST RMP coils
end
if strcmp(Acad.DEVICE,'MAST-U')
dphi = 0.3812;     %MAST-U RMP coils
end
if strcmp(Acad.DEVICE,'AUG')
dphiU = 0.65293; 
dphiL = 0.72955;
dphi = dphiU;     %AUG B-coils
end
end

if n==1
   IC = [1 1 1 -1 -1 -1];
elseif n==2
   IC=[1 -0.5 -0.5 1 -0.5 -0.5];
elseif n==3
   IC = [1 -1 1 -1 1 -1];
end

N   = length(IC);
phi = 2*pi/N*([0:N-1]+0.5);

fA  = n*pi/sin(n*dphi/2)/sum(IC.*exp(-i*n*phi));
if n==2, fA  = n*pi/sin(n*dphi/2)/N;  end %for D3D, n=2 cos-wave form
fAA = abs(fA);

% coil current into MARS-F 
I0  = 1.0e+3;
mu0 = 4.0e-7*pi;
fI  = R0*B0/mu0;
FEEDI = I0/fAA/fI;


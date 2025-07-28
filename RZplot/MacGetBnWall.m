% compute poloidal Fourier harmonics of Bn
% either from MARS-F data BPLASMA: KdataType=1
% or from ERGOS data: KdataType=2

function MacGetBnWall(R,Z,Br,Bz,Bphi)

global Mac
global SDIR

KdataType = 1;

Bnm = 0;  

% find the index for wall surfaces
II = zeros(1,length(Mac.rw));
for j=1:length(II)
    [smin,II(j)] = min(abs(Mac.s-Mac.rw(j)));
end

s_rw=Mac.s(II);

Tw = [];  Rw = [];  Zw = []; Cw = [];
Bwr = []; Bwz = []; Bwphi = [];

if KdataType==1
for j=1:1
    Zmin=min(Z(II(j),:));  Zmax=max(Z(II(j),:)); Z0=0.5*(Zmin+Zmax);
    Rmin=min(R(II(j),:));  Rmax=max(R(II(j),:)); R0=0.5*(Rmin+Rmax);
    R0 = R(1,1); Z0=Z(1,1);
    T = atan2(Z(II(j),:)-Z0,R(II(j),:)-R0);
    [T,JJ] = sort(T);
    Tw = [Tw T];
    Rw = [Rw R(II(j),JJ)];
    Zw = [Zw Z(II(j),JJ)];
    Bwr = [Bwr Br(II(j),JJ)];
    tmp = Bz(II(j),JJ);
    Bwz = [Bwz tmp];
    tmp = Bphi(II(j),JJ); 
    Bwphi = [Bwphi tmp];

    Cw = [Cw Mac.chi(JJ)];
end
end
  
%get Bn, Bt, Bphi in Fourier harmonics in physical theta angle
%using Gauss quadrature for integration
if length(Mac.rw)==1 

%clean Tw,Rw,Zw,Bwr,Bwz,Bwphi
[Tw,JJ] = unique(Tw); Rw = Rw(JJ); Zw = Zw(JJ); Cw = Cw(JJ);
Bwr = Bwr(JJ); Bwz = Bwz(JJ); Bwphi = Bwphi(JJ); 
Tw = [Tw(end)-2*pi Tw]; Rw = [Rw(end) Rw]; Zw=[Zw(end) Zw];
Bwr = [Bwr(end) Bwr];  Bwz = [Bwz(end) Bwz];  Bwphi = [Bwphi(end) Bwphi];
Cw  = [Cw(end)-2*pi Cw];

%get Gauss quadrature points for Fourier decomposition
[z,w] = MacGaussQuad1D(2);
x0 = (Tw(1:end-1)+Tw(2:end))*0.5;
h2 = diff(Tw)*0.5;
xx = z'*h2 + ones(size(z'))*x0;
wh = w'*h2;
xx = xx(:);  wh=wh(:);

kdR = 0;
if kdR > 0 

%Fourier decompose (Rw,Zw) in Tw-angle
expmt = exp(-i*xx*[-Mac.Nm2:Mac.Nm2])/(2*pi);
Rx    = spline(Tw,Rw,xx').*wh';
Zx    = spline(Tw,Zw,xx').*wh';
Rm = Rx*expmt; 
Zm = Zx*expmt;

%compute (dR/dTw,dZ/dTw) via Fourier harmonics
expmt = exp(i*[-Mac.Nm2:Mac.Nm2]'*Tw);
dR = i*([-Mac.Nm2:Mac.Nm2].*Rm)*expmt;
dZ = i*([-Mac.Nm2:Mac.Nm2].*Zm)*expmt;
dR = real(dR);
dZ = real(dZ);

else

%another way to compute dR,dZ
dRT = diff(Rw)./diff(Tw);
x = (Tw(2:end)+Tw(1:end-1))/2;
x = [x(end)-2*pi x x(1)+2*pi];
dRT = [dRT(end) dRT dRT(1)];
dR = pchip(x,dRT,Tw);
dZT = diff(Zw)./diff(Tw);
dZT = [dZT(end) dZT dZT(1)];
dZ = pchip(x,dZT,Tw);

end

%compute physical Bn & Bt from Bwr & Bwz
nA = sqrt(dR.^2+dZ.^2);
Bn = (Bwr.*dZ - Bwz.*dR)./nA;

%Fourier decompose Bn, Bt & Bphi along Cw angle
x0 = (Cw(1:end-1)+Cw(2:end))*0.5;
h2 = diff(Cw)*0.5;
xx = z'*h2 + ones(size(z'))*x0;
wh = w'*h2;
xx = xx(:);  wh=wh(:);

expmt = exp(-i*xx*Mac.Mm')/(2*pi);
Bnx   = spline(Cw,Bn,xx').*wh';
Bnm = Bnx*expmt; 

% save results
% note that this data is for B_n not B^1, cannot be used for ESC procedure
%res = [length(Mac.Mm) II(1); [real(Bnm(:)) imag(Bnm(:))]];
%save INPUT_BNM res -ascii -double
end

% read B1U from B123U data, and if IB>0, 
% convert to Bn in Fourier harmonics in geometrical angle 

function BF = MacPrmBr(B123U,IB,R,dRdchi,dZdchi,T,IIT)

global Mac

NM = round(length(B123U)/6);

if IB==0
   
BF = B123U(1:NM) + B123U(NM+1:2*NM)*i;
BF = BF(:);

else

B1m = B123U(1:NM) + B123U(NM+1:2*NM)*i;
G22  = dRdchi(IB,:).^2 + dZdchi(IB,:).^2;

% compute Bn
tt = Mac.chi';
mm = Mac.Mm';
Fexp = exp(i*tt*mm);
B1 = Fexp * B1m(:);
RIB = R(IB,:);
Bn = B1./RIB(:)./sqrt(G22(:));

% decompose Bn in Fourier space in geometrical angle T
Bn = Bn(IIT);
[z,w] = MacGaussQuad1D(2);
x0 = (T(1:end-1)+T(2:end))*0.5;
h2 = diff(T)*0.5;
xx = z'*h2 + ones(size(z'))*x0;
wh = w'*h2;
xx = xx(:);  wh=wh(:);
expmt = exp(-i*xx*Mac.Mm')/(2*pi);

Bnx   = spline(T,Bn,xx').*wh';
BF = Bnx*expmt; 
BF = BF(:);


end

BF = BF*Mac.BNORM;

 

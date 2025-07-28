% read B2U from B123U data, and if IB>0, 
% convert to Bn in Fourier harmonics 

function BF = MacPrmBt(B123U,IB)

NM = round(length(B123U)/6);

if IB==0
   
BF = B123U(2*NM+1:3*NM) + B123U(3*NM+1:4*NM)*i;

else

BF = B123U(2*NM+1:3*NM) + B123U(3*NM+1:4*NM)*i;

end

BF = BF(:);
BF = BF*Mac.BNORM;


 

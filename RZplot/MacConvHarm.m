% Convert poloidal Fourier harmonics from PEST to EQAC coordinates
% at one rational surface
% Read input angles information from PANGLE_EQAC.IN and PANGLE_PEST.IN

SDIR  = '';
ISURF = 1;   %the number of q=m/n rational surface in PANGLE-files, where to do the conversion
JSURF = 110; %radial mesh point corresponding to the q=m/n rational surface
m     = 2;   %input single m-harmonic in PEST coordinates
M1    =-33;  %output harmonics range [M1,M2] in EQAC coordinates
M2    = 33;   

data_eqac = load([SDIR 'PANGLE_EQAC.IN']);
data_pest = load([SDIR 'PANGLE_PEST.IN']);
 
NRATSURF = data_eqac(1);
NCHIP1   = data_eqac(2);   %NCHI+1
  
eqac = data_eqac(3:end);
pest = data_pest(3:end);

eqac = reshape(eqac,NCHIP1,NRATSURF);
pest = reshape(pest,NCHIP1,NRATSURF);

%note that both eqac and pest are function of uniform geometric angle theta
eqac = eqac(:,ISURF);
pest = pest(:,ISURF);

%assuming a unit single harmonic m in PEST coordinates
%calculate all k-harmonics in EQAC coordinates, in the range of [M1,M2]
b1m = 1.0;
b1k = zeros(M2-M1+1,1);
for k=M1:M2
    f = exp(i*m*pest-i*k*eqac);
    b1k(k-M1+1) = sum((f(1:end-1)+f(2:end)).*diff(eqac));
end
b1k = b1k*b1m/4.0/pi;

%output result
%res = [[M1:M2]' real(b1k) imag(b1k)]
res = [length([M1:M2]) JSURF; real(b1k) imag(b1k)];
save INPUT_BNM.IN res -ascii -double



% New way of generating BNORM01
% minimizing sideband effect at m=M1 and M2

mm = [-9:39]';
M  = length(mm);
dt = 0.5*pi;
t  = linspace(-pi,pi,max(abs(mm))*32+1);
y  = exp(-t.^2/dt^2); 
y  = y - sum(y)*(t(2)-t(1))/3/pi;

%get Gauss quadrature points for Fourier decomposition
[z,w] = MacGaussQuad1D(4);
x0 = (t(1:end-1)+t(2:end))*0.5;
h2 = diff(t)*0.5;
xx = z'*h2 + ones(size(z'))*x0;
wh = w'*h2;
xx = xx(:);  wh=wh(:);

%Fourier decompose y in Tw-angle
expmt = exp(-i*xx*mm')/(2*pi);
yy    = spline(t,y,xx').*wh';
Bnm0  = yy*expmt; 
Bnm0  = Bnm0(:);

for k = 1:M
  Ts  = 2*pi*(k-1)/M;
  Bnm = Bnm0.*exp(i*mm*Ts);
  Bn = [real(Bnm) imag(Bnm)];

  scom = ['save BNORM01n' num2str(k) ' Bn -ascii -double'];
  eval(scom)
end

t = linspace(-pi,pi,201)';
tm = t*mm';
Bn = exp(1i*tm)*Bnm0;
figure(1)
plot(t,real(Bn),'r-',t,imag(Bn),'b--'), hold on,

res = [mm abs(Bnm0)/max(abs(Bnm0))]


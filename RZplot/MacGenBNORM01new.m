% New way of generating BNORM01
% minimizing sideband effect at m=M1 and M2

mm = [-9:39]';
M = length(mm);
%C = 4.3359e-1;
C = 1.0000e+6;
Bnm0 = exp(-(mm-3).^2/max((mm-3).^2)/C);
BGG  = zeros(M,M);
for k = 1:M
  Ts  = 2*pi*(k-1)/M;
  Bnm = Bnm0.*exp(i*mm*Ts);
  BGG(:,k) = Bnm;
  Bn = [real(Bnm) imag(Bnm)];

  scom = ['save BNORM01n' num2str(k) ' Bn -ascii -double'];
  eval(scom)
end

disp(['Condition number = ', num2str(cond(BGG))])  %=1.2561e+01

t = linspace(-pi,pi,1001)';
tm = t*mm';
Bn = exp(1i*tm)*Bnm0;
figure(1)
plot(t,real(Bn),'r-',t,imag(Bn),'b--'), hold on,

res = [mm abs(Bnm0)/max(abs(Bnm0))]


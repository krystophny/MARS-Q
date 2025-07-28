function [JM1C,JM2C,JM3C] = MacGetJM123(CHIG,WG,J1G,J2G,J3G)

global Mac

JM1C = zeros(size(Mac.Mm));
JM2C = zeros(size(Mac.Mm));
JM3C = zeros(size(Mac.Mm));

for m=1:length(Mac.Mm)
  expchiG = exp(-i*Mac.Mm(m)*CHIG);
  JM1C(m) = sum(WG.*sum(J1G.*expchiG,2));
  JM2C(m) = sum(WG.*sum(J2G.*expchiG,2));
  JM3C(m) = sum(WG.*sum(J3G.*expchiG,2));
end

h = Mac.chi(2) - Mac.chi(1);
JM1C = JM1C*h/4/pi;
JM2C = JM2C*h/4/pi;
JM3C = JM3C*h/4/pi;

%recompute JM3 using DivJ=0
JM3Cn = -Mac.Mm.*JM2C/Mac.n;
JM3Ce = abs(JM3C-JM3Cn)./abs(JM3C);

%show results
%JM12C = [JM1C JM2C]
%JM33C = [JM3C JM3Cn JM3Ce]

%direct DivJ=0 check
J1 = 0; J2 = 0; J3 = 0; DivJ = 0;
for m=1:Mac.Nm1
  mm = Mac.Mm(m);
  em = exp(i*mm*Mac.chi);
  J1 = J1 + JM1C(m)*em;
  J2 = J2 + JM2C(m)*em;
  J3 = J3 + JM3C(m)*em;
  DivJ = DivJ + i*(mm*JM2C(m)+Mac.n*JM3C(m))*em;
end

%FEEDJ0 = real(sum(JM2C))

if Mac.plot_shdiv > 0
  figure(Mac.plot_shdiv)
  subplot(2,2,1), plot(Mac.chi,real(J1),'r-','LineWidth',2), hold on,
                  plot(Mac.chi,imag(J1),'r--','LineWidth',2), hold on,
                  ylabel('J^1','FontSize',16)
  subplot(2,2,2), plot(Mac.chi,real(J2),'r-','LineWidth',2), hold on,
                  plot(Mac.chi,imag(J2),'r--','LineWidth',2), hold on,
                  ylabel('J^2','FontSize',16)
  subplot(2,2,3), plot(Mac.chi,real(J3),'r-','LineWidth',2), hold on,
                  plot(Mac.chi,imag(J3),'r--','LineWidth',2), hold on,
                  xlabel('\chi','FontSize',16)
                  ylabel('J^3','FontSize',16)
  subplot(2,2,4), plot(Mac.chi,real(DivJ),'r-','LineWidth',2), hold on,
                  plot(Mac.chi,imag(DivJ),'r--','LineWidth',2), hold on,
                  xlabel('\chi','FontSize',16)
                  ylabel('DivJ','FontSize',16)
end

%save data
tmp = [Mac.Mm real(JM2C) imag(JM2C)];
fid = fopen('FEEDJ_CARIDDI','w');
fprintf(fid,'%4i\n',Mac.Nm1);
fprintf(fid,'%3i  %23.16e  %23.16e\n',tmp');
fclose(fid);


function W = MacGetNTVW

global Mac

W  = zeros(Mac.Ns1,Mac.Nm1);
kk = linspace(0,1,101).^2;
[K,E] = ellipke(kk);
KE = ones(size(Mac.q))*(E-(1-kk).*K);
KE(:,1) = KE(:,2);  KE(:,end) = 1;

for mm=1:length(Mac.Mm)
  m = Mac.Mm(mm);
  Fmc = zeros(Mac.Ns1,length(kk));
  Fms = Fmc;
  for j=2:length(kk) 
    k2 = kk(j);
    tt = linspace(0,2*asin(sqrt(k2)),101);
    FF = ones(size(Mac.q))*sqrt(abs(k2-sin(tt/2).^2));
    Fc = cos((m+Mac.n*Mac.q)*tt).*FF;
    Fs = sin((m+Mac.n*Mac.q)*tt).*FF;
    Fmc(:,j) = (sum(Fc,2)-Fc(:,1)*0.5-Fc(:,end)*0.5)*tt(2)*2;
    Fms(:,j) = (sum(Fs,2)-Fs(:,1)*0.5-Fs(:,end)*0.5)*tt(2)*2;
  end

  F = (Fmc.^2+Fms.^2)./KE;
  G = (F(:,1:end-1)+F(:,2:end))*0.5.*(ones(size(Mac.q))*diff(kk));
  W(:,mm) = sum(G,2);
end

if Mac.plot_NTVW
  figure(Mac.plot_NTVW)
  subplot(2,1,1), plot(Mac.s(1:Mac.Ns1),W), hold on,
                  xlabel('s','FontSize',16),
                  ylabel('W_{nm}','FontSize',16),
  subplot(2,1,2), yy = abs(W).^2; 
                  plot(Mac.Mm,sqrt(sum(yy,1)/size(yy,1))), hold on,
                  ylabel('||W_{nm}||_{L_2}','FontSize',14),
                  xlabel('m','FontSize',14),
end


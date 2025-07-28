function Torq = MacGetJXBTorq(Rmid,JM1,JM2,BM1,BM2);

global Mac

N = Mac.Ns1;
tmp = conj(JM1(1:N,:)).*BM2(1:N,:) - conj(JM2(1:N,:)).*BM1(1:N,:);
Torq = Mac.R0EXP^2*2*pi^2*real(sum(tmp,2));

if Mac.plot_JXB > 0
  figure(Mac.plot_JXB)
  plot(Rmid,Torq,'b--','LineWidth',2), hold on,
  xlabel('R [m]','FontSize',16),
  ylabel('torque density [N]','FontSize',16),
end

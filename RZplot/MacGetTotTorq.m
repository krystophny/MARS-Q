function Torq = MacGetTotTorq(Rmid,TorqNTV,TorqJXB);

global Mac

Torq = TorqNTV + TorqJXB;

if Mac.plot_Torq > 0
  figure(Mac.plot_Torq)
  plot(Rmid,Torq,'r-','LineWidth',2), hold on,
  xlabel('R [m]','FontSize',16),
  ylabel('torque density [N]','FontSize',16),
end

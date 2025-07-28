function Sv = MacGetNTVSv(Rmid,jacobian)

global Mac

jacob = jacobian(1:Mac.Ns1,:);
Sv = 2*pi*Mac.R0EXP^2*(sum(jacob,2)-jacob(:,1)/2-jacob(:,end)/2)*(Mac.chi(2)-Mac.chi(1));

if Mac.plot_NTVSv
  figure(Mac.plot_NTVSv)
  plot(Rmid,Sv,'b-','LineWidth',2), hold on,
  xlabel('R [m]','FontSize',16),
  ylabel('innertial shell volume [m^2]','FontSize',14),
end

data = [Rmid Sv];
save dataSv data -ascii -double

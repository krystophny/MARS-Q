function MacPlotBeq(jacobian,dpsids,dRdchi,dZdchi)

global Mac

ssplot  = [0.5 0.7 0.8];
sscol   = ['rbk'];

hf=figure(Mac.plot_BEQ);

for k=1:length(ssplot)
splot = ssplot(k);
[Y,II] = min(abs(Mac.s-splot));
BR = dpsids(II)./jacobian(II,:).*dRdchi(II,:);
plot(Mac.chi*180/pi,BR*Mac.B0EXP,'-','Color',sscol(k),'LineWidth',2), hold on,
end

xlabel('equal-arc pol. angle [deg]','FontSize',18)
ylabel('B_R [T]','FontSize',18)
ha=get(hf,'CurrentAxes');
set(ha,'FontSize',18)
legend(num2str(ssplot(1)),num2str(ssplot(2)),num2str(ssplot(3)))






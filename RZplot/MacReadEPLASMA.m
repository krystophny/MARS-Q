function MacReadEPLASMA(filename)

global Mac

EPLASMA = load(filename);

s   = EPLASMA(:,1);
dWI = EPLASMA(:,5) + EPLASMA(:,6)*i;
dWF = EPLASMA(:,7) + EPLASMA(:,8)*i;
dWK = EPLASMA(:,9) + EPLASMA(:,10)*i;

RES_DW = [sum(dWI(4:end)) sum(dWF(4:end)) sum(dWK(4:end))]
fac_wk = -2.8071;

%central point
dWI(1) = 0;
dWF(1) = 0;
dWK(1) = 0;

if Mac.plot_ER > 0
   figure(Mac.plot_ER)
   plot(s,real(dWI),'b-',s,real(dWF),'b--',s,real(dWK)*fac_wk,'r-.','LineWidth',2), hold on,
   ylabel('dW/d\psi_p','FontSize',16,'FontWeight','Bold')
   xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha = get(Mac.plot_ER,'CurrentAxes');
   set(ha,'FontSize',14,'FontWeight','Bold')
   axis tight
   legend('fluid','para.kin.','perp.kin.')
end

if Mac.plot_EI > 0
   figure(Mac.plot_EI)
   plot(s,imag(dWI),'b-',s,imag(dWF),'b--',s,imag(dWK),'r-.','LineWidth',2), hold on,
   ylabel('dW/d\psi_p','FontSize',16,'FontWeight','Bold')
   xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha = get(Mac.plot_EI,'CurrentAxes');
   set(ha,'FontSize',14,'FontWeight','Bold')
   legend('fluid','para.kin.','perp.kin.')
end

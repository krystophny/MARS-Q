function [TC, BrC, BzC, BpC] = MacReadBCARIDDI(filename)

global Mac

command = ['!cp ' filename ' MacDataBs'];
eval(command), load MacDataBs,

RC = MacDataBs(:,1);
ZC = MacDataBs(:,2);

fac = 1.0/4.0e-7/pi;
BrC = fac*(MacDataBs(:,3)+MacDataBs(:,9) + (MacDataBs(:,4)+MacDataBs(:,10))*i);
BzC = fac*(MacDataBs(:,5)+MacDataBs(:,11) + (MacDataBs(:,6)+MacDataBs(:,12))*i);
BpC =-fac*(MacDataBs(:,7)+MacDataBs(:,13) + (MacDataBs(:,8)+MacDataBs(:,14))*i);
TC = atan2(ZC,RC-Mac.R0EXP);

[TC,II,JJ] = unique(TC);
BrC = BrC(II); BzC = BzC(II); BpC = BpC(II); 

if Mac.plot_BC > 0
  figure(Mac.plot_BC)
  subplot(3,2,1), plot(TC,real(BrC),'b-','LineWidth',2), hold on
                  ylabel('Re(Br)','FontSize',14) 
  subplot(3,2,2), plot(TC,imag(BrC),'b-','LineWidth',2), hold on 
                  ylabel('Im(Br)','FontSize',14) 
  subplot(3,2,3), plot(TC,real(BzC),'b-','LineWidth',2), hold on 
                  ylabel('Re(Bz)','FontSize',14) 
  subplot(3,2,4), plot(TC,imag(BzC),'b-','LineWidth',2), hold on 
                  ylabel('Im(Bz)','FontSize',14) 
  subplot(3,2,5), plot(TC,real(BpC),'b-','LineWidth',2), hold on 
                  ylabel('Re(B_\phi)','FontSize',14) 
                  xlabel('\theta','FontSize',14) 
  subplot(3,2,6), plot(TC,imag(BpC),'b-','LineWidth',2), hold on 
                  ylabel('Im(B_\phi)','FontSize',14) 
                  xlabel('\theta','FontSize',14) 
end

res = [TC real(BrC) imag(BrC) real(BzC) imag(BzC) real(BpC) imag(BpC)];
save Bfield_CARIDDI res -ascii 

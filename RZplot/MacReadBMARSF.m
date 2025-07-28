function [TC, BrC, BzC, BpC] = MacReadBMARSF(filename)

global Mac

command = ['!cp ' filename ' MacDataBs'];
eval(command), load MacDataBs,

RC = MacDataBs(:,1);
ZC = MacDataBs(:,2);

BrC = MacDataBs(:,3) + MacDataBs(:,4)*i;
BzC = MacDataBs(:,5) + MacDataBs(:,6)*i;
BpC = MacDataBs(:,7) + MacDataBs(:,8)*i;

TC = atan2(ZC,RC-Mac.R0EXP);

[TC,II,JJ] = unique(TC);
BrC = BrC(II); BzC = BzC(II); BpC = BpC(II); 

if Mac.plot_BC == 10
  figure(Mac.plot_BC)
  subplot(3,2,1), plot(TC,real(BrC),'b-','LineWidth',2), hold on 
  subplot(3,2,2), plot(TC,imag(BrC),'b-','LineWidth',2), hold on 
  subplot(3,2,3), plot(TC,real(BzC),'b-','LineWidth',2), hold on 
  subplot(3,2,4), plot(TC,imag(BzC),'b-','LineWidth',2), hold on 
  subplot(3,2,5), plot(TC,real(BpC),'b-','LineWidth',2), hold on 
  subplot(3,2,6), plot(TC,imag(BpC),'b-','LineWidth',2), hold on 
end

res = [TC real(BrC) imag(BrC) real(BzC) imag(BzC) real(BpC) imag(BpC)];
save Bfield_MARSF res -ascii 

function [TC, JrC, JzC, JpC] = MacReadJCARIDDI(filename)

global Mac

command = ['!cp ' filename ' MacDataJs'];
eval(command), load MacDataJs,

RC = MacDataJs(:,1);
ZC = MacDataJs(:,2);

fac = 4.0e-7*pi*Mac.R0EXP/Mac.B0EXP/Mac.HWEXP/Mac.ALPHA;
JrC = fac*(MacDataJs(:,3)+MacDataJs(:,4)*i);
JzC = fac*(MacDataJs(:,5)+MacDataJs(:,6)*i);
JpC =-fac*(MacDataJs(:,7)+MacDataJs(:,8)*i);
TC = atan2(ZC,RC-Mac.R0EXP);

[TC,II,JJ] = unique(TC);
JrC = JrC(II); JzC = JzC(II); JpC = JpC(II); 

if Mac.plot_FEEDJ > 0
  figure(Mac.plot_FEEDJ)
  subplot(3,2,1), plot(TC,real(JrC),'b-','LineWidth',2), hold on
                  ylabel('Re(Jr)','FontSize',14) 
  subplot(3,2,2), plot(TC,imag(JrC),'b-','LineWidth',2), hold on 
                  ylabel('Im(Jr)','FontSize',14) 
  subplot(3,2,3), plot(TC,real(JzC),'b-','LineWidth',2), hold on 
                  ylabel('Re(Jz)','FontSize',14) 
  subplot(3,2,4), plot(TC,imag(JzC),'b-','LineWidth',2), hold on 
                  ylabel('Im(Jz)','FontSize',14) 
  subplot(3,2,5), plot(TC,real(JpC),'b-','LineWidth',2), hold on 
                  ylabel('Re(J_\phi)','FontSize',14) 
                  xlabel('\theta','FontSize',14) 
  subplot(3,2,6), plot(TC,imag(JpC),'b-','LineWidth',2), hold on 
                  ylabel('Im(J_\phi)','FontSize',14) 
                  xlabel('\theta','FontSize',14) 
end


% plot J at the shell
if Mac.plot_shell>0
   figure(Mac.plot_shell)
   h = 6.0/Mac.Ns1;
   Jrr = real(JrC); Jzr = real(JzC);
   Jt = sqrt(Jrr.^2 + Jzr.^2);  
   %Jt = max(max(Jt));
   [II,JJ]=find(Jt==0.0);  Jt(II,JJ) = 1.0;
   R1 = RC + h*Jrr./Jt;
   Z1 = ZC + h*Jzr./Jt;

   R2 = [RC(:) R1(:)]';
   Z2 = [ZC(:) Z1(:)]';

   plot(RC,ZC,'c-'), hold on,
   plot(R2,Z2,'b-'), hold on,
   axis equal
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
end



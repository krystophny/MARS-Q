function [J1,J2,J3,JRE] = MacGetJ123(JM1,JM2,JM3,JMRE,R,Z)

global Mac

expmchi = exp(Mac.Mm*Mac.chi*i);

J1 = JM1*expmchi;
J2 = JM2*expmchi;
J3 = JM3*expmchi;
JRE = JMRE*expmchi;


%plot current on the shell surface
if Mac.plot_Js > 0
  fac = Mac.ShellCurrNorm;
  [smin,II] = min(abs(Mac.s-Mac.rs(1)));
  DivJ = i*((Mac.Mm)'.*JM2(II,:)+Mac.n*JM3(II,:))*expmchi;
  AmpJ20 = fac*real(sum(JM2(II,:)))
  figure(Mac.plot_Js)
% subplot(2,2,1), plot(Mac.chi,real(fac*J1(II,:)),'b-','LineWidth',2), hold on,
%                 plot(Mac.chi,imag(fac*J1(II,:)),'b--','LineWidth',2), hold on,
%                 ylabel('J^1','FontSize',16)
  subplot(2,1,1), plot(Mac.chi,real(fac*J2(II,:)),'b-','LineWidth',2), hold on,
                  plot(Mac.chi,imag(fac*J2(II,:)),'b--','LineWidth',2), hold on,
                  ylabel('J^2','FontSize',16)
  subplot(2,1,2), plot(Mac.chi,real(fac*J3(II,:)),'b-','LineWidth',2), hold on,
                  plot(Mac.chi,imag(fac*J3(II,:)),'b--','LineWidth',2), hold on,
                  xlabel('\chi','FontSize',16)
                  ylabel('J^3','FontSize',16)
% subplot(2,2,4), plot(Mac.chi,real(fac*DivJ),'b-','LineWidth',2), hold on,
%                 plot(Mac.chi,imag(fac*DivJ),'b--','LineWidth',2), hold on,
%                 xlabel('\chi','FontSize',16)
%                 ylabel('DivJ','FontSize',16)
end

if Mac.plot_J2J3 > 0
   hf=figure(10*Mac.plot_J2J3+0);
   [smin,II] = min(abs(Mac.s-Mac.rw(1)));
   Rw = R(II,:);
   Zw = Z(II,:);
   Tw = atan2(Zw-Z(1,1),Rw-R(1,1))*180/pi;
   [Tw,JJ]=sort(Tw);

   hs=subplot(2,2,1);
   plot(Tw,real(J2(II,JJ)),'b-','LineWidth',3), hold on,
   ylabel('Re(J^2)','FontSize',16,'FontWeight','Bold')
   set(hs,'FontSize',14,'FontWeight','Bold')
   
   hs=subplot(2,2,2);
   plot(Tw,imag(J2(II,JJ)),'b-','LineWidth',3), hold on,
   ylabel('Im(J^2)','FontSize',16,'FontWeight','Bold')
   set(hs,'FontSize',14,'FontWeight','Bold')
   
   hs=subplot(2,2,3);
   plot(Tw,real(J3(II,JJ)),'b-','LineWidth',3), hold on,
   xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
   ylabel('Re(J^3)','FontSize',16,'FontWeight','Bold')
   set(hs,'FontSize',14,'FontWeight','Bold')
   
   hs=subplot(2,2,4);
   plot(Tw,imag(J3(II,JJ)),'b-','LineWidth',3), hold on,
   xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
   ylabel('Im(J^3)','FontSize',16,'FontWeight','Bold')
   set(hs,'FontSize',14,'FontWeight','Bold')
end   

  


function [JM1,JM2,JM3,JMRE] = MacReadJPLASMA(filename,R,Z,jacobian)

global Mac

JPLASMA = load(filename);

Mac.Nm1 = JPLASMA(1,1);

if Mac.Ns ~= JPLASMA(1,2),
    disp('Number of radial points is different for equilibrium and stability!')
   Ns = JPLASMA(1,2);
   Mac.Ns = Ns;
   Mac.s  = Mac.s(1:Ns);
else 
   Ns = JPLASMA(1,2);
end

Mac.n = round(JPLASMA(1,3));
Mac.Mm = round(JPLASMA(2:Mac.Nm1+1,1));

if length(Mac.mm_plot)>0 
   Mac.mm_true = Mac.mm_plot-Mac.Mm(1)+1;
else
   Mac.mm_true = Mac.Mm-Mac.Mm(1)+1;
end

JM1 = JPLASMA(Mac.Nm1+2:end,1) + JPLASMA(Mac.Nm1+2:end,2)*i;
JM2 = JPLASMA(Mac.Nm1+2:end,3) + JPLASMA(Mac.Nm1+2:end,4)*i;
JM3 = JPLASMA(Mac.Nm1+2:end,5) + JPLASMA(Mac.Nm1+2:end,6)*i;
if size(JPLASMA,2)==8
   JMRE = JPLASMA(Mac.Nm1+2:end,7) + JPLASMA(Mac.Nm1+2:end,8)*i;
elseif size(JPLASMA,2)==14
   JM1 = JPLASMA(Mac.Nm1+2:end,7) + JPLASMA(Mac.Nm1+2:end,8)*i;
   JM2 = JPLASMA(Mac.Nm1+2:end,9) + JPLASMA(Mac.Nm1+2:end,10)*i;
   JM3 = JPLASMA(Mac.Nm1+2:end,11) + JPLASMA(Mac.Nm1+2:end,12)*i;

   JMRE = JPLASMA(Mac.Nm1+2:end,13) + JPLASMA(Mac.Nm1+2:end,14)*i;
else
   JMRE = JM3;
end

JM1 = reshape(JM1,Mac.Ns,Mac.Nm1);
JM2 = reshape(JM2,Mac.Ns,Mac.Nm1);
JM3 = reshape(JM3,Mac.Ns,Mac.Nm1);
JMRE = reshape(JMRE,Mac.Ns,Mac.Nm1);


JM1 = JM1*Mac.JNORM; 
JM2 = JM2*Mac.JNORM; 
JM3 = JM3*Mac.JNORM; 
JMRE = JMRE*Mac.JNORM; 

% save JM2 to CURHARMO.OUT
if 1==0
II = 225;
res = [Mac.Nm1 Mac.Nm1 Mac.Nm1; Mac.Mm real(JM2(II,:))' imag(JM2(II,:))'];
save CURHARMO.OUT res -ascii -double
end

%m0 = -Mac.Mm(1) + 1;
%JM2(:,m0) = 0.0;

% note that JM1 is defined at half-points, recompute at integer-points
if Mac.spline_J1
  x = (Mac.s(1:Mac.Ns-1) + Mac.s(2:Mac.Ns))*0.5;
  JM1new = JM1;
  JM1new(2:end-1,:) = transpose(pchip(x',transpose(JM1(1:end-1,:)),Mac.s(2:Mac.Ns-1)'));
  JM1new(1,:) = 0;  JM1new(end,:) = 0;
  JM1 = JM1new;
  clear JM1new
end

% patch first three points
JM1(1:3,:) = 0;
JM2(1:3,:) = 0;
JM3(1:3,:) = 0;

%remove surface currents, should be avoided for MacSolTest.m 
if 1==0
  II = find(Mac.s>0.9985 & Mac.s<=1.0); I1 = II(1)-1;
  for k=1:size(JM1,2)
  JM1(II,k) = JM1(I1,k);
  JM2(II,k) = JM2(I1,k);
  JM3(II,k) = JM3(I1,k);
  end
end

%  get phase of JM2 for m=2 harmonic at rationa surface Iratsurf(1)
%  and modify phase to remove real part
if 1==0
   m0 = 2;
   I0 = Mac.Iratsurf(1);
   m1 = m0 - Mac.Mm(1) + 1;
   p0 =-angle(JM2(I0,m1))+pi/2;
   f0 = exp(i*p0);
   JM1 = JM1*f0;
   JM2 = JM2*f0;
   JM3 = JM3*f0;
end
     
if Mac.plot_JM > 0
   figure(10*Mac.plot_JM+0)
   SS = '-';
   subplot(4,2,1), plot(Mac.s(1:Mac.Ns),real(JM1),SS), hold on,
                   ylabel('Re(J^1_m)','FontSize',14)
   subplot(4,2,2), plot(Mac.s(1:Mac.Ns),imag(JM1),SS), hold on,
                   xlabel('s','FontSize',14), ylabel('Im(J^1_m)','FontSize',14)

   subplot(4,2,3), plot(Mac.s(1:Mac.Ns),real(JM2),SS), hold on,
                   ylabel('Re(J^2_m)','FontSize',14)
   subplot(4,2,4), plot(Mac.s(1:Mac.Ns),imag(JM2),SS), hold on,
                   xlabel('s','FontSize',14), ylabel('Im(J^2_m)','FontSize',14)

   subplot(4,2,5), plot(Mac.s(1:Mac.Ns),real(JM3),SS), hold on,
                   ylabel('Re(J^3_m)','FontSize',14)
   subplot(4,2,6), plot(Mac.s(1:Mac.Ns),imag(JM3),SS), hold on,
                   xlabel('s','FontSize',14), ylabel('Im(J^3_m)','FontSize',14)

   subplot(4,2,7), plot(Mac.s(1:Mac.Ns),real(JMRE),SS), hold on,
                   ylabel('Re(J_{RE})','FontSize',14)
   subplot(4,2,8), plot(Mac.s(1:Mac.Ns),imag(JMRE),SS), hold on,
                   xlabel('s','FontSize',14), ylabel('Im(J_{RE})','FontSize',14)

   sleft = 0;
   mleft = 0;
   plot_JM_EQ = Mac.plot_JM_EQ; 
                   %=1: plot equilibrium current density
   plot_JM_KD = Mac.plot_JM_KD; 
                   %=1: 1-D plot along s for Fourier harmonics; 
                   %=2: 1-D plot along m for all s
                   %=3: 2-D plot in (m,s)-space
                   %=4: 2-D plot in (chi,s)-space
		   %=5: 2-D plot in (R,Z)-space
   ssm = (Mac.s(1:Mac.Ns1-1)+Mac.s(2:Mac.Ns1))/2;

   if plot_JM_KD==1
   hf=figure(10*Mac.plot_JM+1);
   plot(ssm,real(JM1(1:Mac.Ns1-1,:)),SS,'LineWidth',2), hold on,
   xlabel('s','FontSize',18)
   ylabel('Re(J^1_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([sleft 1 a(3) a(4)]);
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end

   hf=figure(10*Mac.plot_JM+2);
   plot(ssm,real(JM1(1:Mac.Ns1-1,:)),SS,'LineWidth',2), hold on,
   xlabel('s','FontSize',18)
   ylabel('Re(J^1_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([sleft 1 a(3) a(4)]);
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end

   hf=figure(10*Mac.plot_JM+3);
   plot(Mac.s(1:Mac.Ns1),real(JM2(1:Mac.Ns1,:)),SS,'LineWidth',2), hold on,

   %plot equilibrium current density
   if plot_JM_EQ==1
   JEQ2 =-diff(Mac.TEQ)./diff(Mac.s(1:Mac.Ns1));
   sx   = (Mac.s(1:Mac.Ns1-1)+Mac.s(2:Mac.Ns1))/2;
   plot(sx,-JEQ2,'b--','LineWidth',2), hold on,
   end

   xlabel('s','FontSize',18)
   ylabel('Re(J^2_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([sleft 1 a(3) a(4)]);
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end

   hf=figure(10*Mac.plot_JM+4);
   plot(Mac.s(1:Mac.Ns1),imag(JM2(1:Mac.Ns1,:)),SS,'LineWidth',2), hold on,
   xlabel('s','FontSize',18)
   ylabel('Im(J^2_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([sleft 1 a(3) a(4)]);
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end

   hf=figure(10*Mac.plot_JM+5);
   plot(Mac.s(1:Mac.Ns1),real(JM3(1:Mac.Ns1,:)),SS,'LineWidth',2), hold on,
   xlabel('s','FontSize',18)
   ylabel('Re(J^3_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([sleft 1 a(3) a(4)]);
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end

   hf=figure(10*Mac.plot_JM+6);
   plot(Mac.s(1:Mac.Ns1),imag(JM3(1:Mac.Ns1,:)),SS,'LineWidth',2), hold on,
   xlabel('s','FontSize',18)
   ylabel('Im(J^3_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([sleft 1 a(3) a(4)]);
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end

   mm = Mac.mm_true;
   hf=figure(10*Mac.plot_JM+7);
   plot(Mac.s(1:Mac.Ns1),real(JMRE(1:Mac.Ns1,mm)),SS,'LineWidth',2), hold on,
   xlabel('s','FontSize',18)
   ylabel('Re(J_{RE,mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([sleft 1 a(3) a(4)]);
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end

   hf=figure(10*Mac.plot_JM+8);
   plot(Mac.s(1:Mac.Ns1),imag(JMRE(1:Mac.Ns1,mm)),SS,'LineWidth',2), hold on,
   xlabel('s','FontSize',18)
   ylabel('Im(J_{RE,mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([sleft 1 a(3) a(4)]);
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   end

   if plot_JM_KD==2
   hf=figure(10*Mac.plot_JM+1);
   plot(Mac.Mm,real(JM1(1:Mac.Ns1-1,:)),SS,'LineWidth',2), hold on,
   xlabel('m','FontSize',18)
   ylabel('Re(J^1_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([mleft a(2) a(3) a(4)]);

   hf=figure(10*Mac.plot_JM+2);
   plot(Mac.Mm,imag(JM1(1:Mac.Ns1-1,:)),SS,'LineWidth',2), hold on,
   xlabel('m','FontSize',18)
   ylabel('Im(J^1_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([mleft a(2) a(3) a(4)]);

   hf=figure(10*Mac.plot_JM+3);
   plot(Mac.Mm,real(JM2(1:Mac.Ns1,:)),SS,'LineWidth',2), hold on,
   xlabel('m','FontSize',18)
   ylabel('Re(J^2_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([mleft a(2) a(3) a(4)]);

   hf=figure(10*Mac.plot_JM+4);
   plot(Mac.Mm,imag(JM2(1:Mac.Ns1,:)),SS,'LineWidth',2), hold on,
   xlabel('m','FontSize',18)
   ylabel('Im(J^2_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([mleft a(2) a(3) a(4)]);

   hf=figure(10*Mac.plot_JM+5);
   plot(Mac.Mm,real(JM3(1:Mac.Ns1,:)),SS,'LineWidth',2), hold on,
   xlabel('m','FontSize',18)
   ylabel('Re(J^3_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([mleft a(2) a(3) a(4)]);

   hf=figure(10*Mac.plot_JM+6);
   plot(Mac.Mm,imag(JM3(1:Mac.Ns1,:)),SS,'LineWidth',2), hold on,
   xlabel('m','FontSize',18)
   ylabel('Im(J^3_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   a=axis; axis([mleft a(2) a(3) a(4)]);
   end

   if plot_JM_KD==3
   hf=figure(10*Mac.plot_JM+1);
   pcolor(Mac.Mm,ssm.^2,real(JM1(1:Mac.Ns1-1,:))), colorbar, shading interp,
   xlabel('m','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Re(J^1_{mn})','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([mleft a(2) sleft 1]); %colormap(hot)
   print(10*Mac.plot_JM+1,'-djpeg',[Mac.SDIR '../Resu3DEQ/J1Mre_' Mac.case]); 

   hf=figure(10*Mac.plot_JM+2);
   pcolor(Mac.Mm,ssm.^2,imag(JM1(1:Mac.Ns1-1,:))), colorbar, shading interp,
   xlabel('m','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Im(J^1_{mn})','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([mleft a(2) sleft 1]); %colormap(hot)
   print(10*Mac.plot_JM+2,'-djpeg',[Mac.SDIR '../Resu3DEQ/J1Mim_' Mac.case]); 

   hf=figure(10*Mac.plot_JM+3);
   pcolor(Mac.Mm,Mac.s(1:Mac.Ns1).^2,real(JM2(1:Mac.Ns1,:))), colorbar, shading interp,
   xlabel('m','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Re(J^2_{mn})','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([mleft a(2) sleft 1]); %colormap(hot)
   print(10*Mac.plot_JM+3,'-djpeg',[Mac.SDIR '../Resu3DEQ/J2Mre_' Mac.case]); 

   hf=figure(10*Mac.plot_JM+4);
   pcolor(Mac.Mm,Mac.s(1:Mac.Ns1).^2,imag(JM2(1:Mac.Ns1,:))), colorbar, shading interp,
   xlabel('m','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Im(J^2_{mn})','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([mleft a(2) sleft 1]); %colormap(hot)
   print(10*Mac.plot_JM+4,'-djpeg',[Mac.SDIR '../Resu3DEQ/J2Mim_' Mac.case]); 

   hf=figure(10*Mac.plot_JM+5);
   pcolor(Mac.Mm,Mac.s(1:Mac.Ns1).^2,real(JM3(1:Mac.Ns1,:))), colorbar, shading interp,
   xlabel('m','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Re(J^3_{mn})','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([mleft a(2) sleft 1]); %colormap(hot)
   print(10*Mac.plot_JM+5,'-djpeg',[Mac.SDIR '../Resu3DEQ/J3Mre_' Mac.case]); 

   hf=figure(10*Mac.plot_JM+6);
   pcolor(Mac.Mm,Mac.s(1:Mac.Ns1).^2,imag(JM3(1:Mac.Ns1,:))), colorbar, shading interp,
   xlabel('m','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Im(J^3_{mn})','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([mleft a(2) sleft 1]); %colormap(hot)
   print(10*Mac.plot_JM+6,'-djpeg',[Mac.SDIR '../Resu3DEQ/J3Mim_' Mac.case]); 
   end

   if plot_JM_KD==4 | plot_JM_KD==5
      expmchi = exp(Mac.Mm*Mac.chi*i);
      J1 = JM1*expmchi; 
      J2 = JM2*expmchi;
      J3 = JM3*expmchi;
      JR = JMRE*expmchi;
      NRES = Mac.NRES;
      for k=1:NRES
          JR(k,:) = JR(NRES+1,:);
      end

      %calculate surface averaged JRE
      N   = Mac.Ns1;
      II  = 2:N;
      JRS = Mac.s(1:N)*0;
      JRS(2:N) = sum(JR(II,:).*jacobian(II,:),2)./sum(jacobian(II,:),2);
      JRS(1)   = JRS(2);
   end

   if plot_JM_KD==4
   hf=figure(10*Mac.plot_JM+1);
   pcolor(Mac.chi*180/pi,ssm.^2,real(J1(1:Mac.Ns1-1,:))), colorbar, shading interp,
   xlabel('\chi [degree]','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Re(J^1)','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)
   %print(10*Mac.plot_JM+1,'-djpeg',[Mac.SDIR '../Resu3DEQ/J1re_' Mac.case]); 

   hf=figure(10*Mac.plot_JM+2);
   pcolor(Mac.chi*180/pi,ssm.^2,imag(J1(1:Mac.Ns1-1,:))), colorbar, shading interp,
   xlabel('\chi [degree]','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Im(J^1)','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)
   %print(10*Mac.plot_JM+2,'-djpeg',[Mac.SDIR '../Resu3DEQ/J1im_' Mac.case]); 

   hf=figure(10*Mac.plot_JM+3);
   pcolor(Mac.chi*180/pi,Mac.s(1:Mac.Ns1).^2,real(J2(1:Mac.Ns1,:))), colorbar, shading interp,
   xlabel('\chi [degree]','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Re(J^2)','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)
   %print(10*Mac.plot_JM+3,'-djpeg',[Mac.SDIR '../Resu3DEQ/J2re_' Mac.case]); 

   hf=figure(10*Mac.plot_JM+4);
   pcolor(Mac.chi*180/pi,Mac.s(1:Mac.Ns1).^2,imag(J2(1:Mac.Ns1,:))), colorbar, shading interp,
   xlabel('\chi [degree]','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Im(J^2)','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)
   %print(10*Mac.plot_JM+4,'-djpeg',[Mac.SDIR '../Resu3DEQ/J2im_' Mac.case]); 

   hf=figure(10*Mac.plot_JM+5);
   pcolor(Mac.chi*180/pi,Mac.s(1:Mac.Ns1).^2,real(J3(1:Mac.Ns1,:))), colorbar, shading interp,
   xlabel('\chi [degree]','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Re(J^3)','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)
   %print(10*Mac.plot_JM+5,'-djpeg',[Mac.SDIR '../Resu3DEQ/J3re_' Mac.case]); 

   hf=figure(10*Mac.plot_JM+6);
   pcolor(Mac.chi*180/pi,Mac.s(1:Mac.Ns1).^2,imag(J3(1:Mac.Ns1,:))), colorbar, shading interp,
   xlabel('\chi [degree]','FontSize',16)
   ylabel('\psi_p','FontSize',16)
   title('Im(J^3)','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)
   %print(10*Mac.plot_JM+6,'-djpeg',[Mac.SDIR '../Resu3DEQ/J3im_' Mac.case]); 
   end

   if plot_JM_KD==5
   N     = Mac.Ns1;
   fac_J = Mac.B0EXP/Mac.R0EXP/(4e-7*pi)/1e+6;

   hf=figure(10*Mac.plot_JM+1);
   pcolor(R(1:N,:)*Mac.R0EXP,Z(1:N,:)*Mac.R0EXP,real(J1(1:N,:))*fac_J), colorbar, shading interp,
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
   title('J^1 [MA/m^2]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   axis equal, colormap(jet)
   %a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)

   hf=figure(10*Mac.plot_JM+2);
   pcolor(R(1:N,:)*Mac.R0EXP,Z(1:N,:)*Mac.R0EXP,imag(J1(1:N,:))*fac_J), colorbar, shading interp,
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
   title('Im(J^1) [MA/m^2]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   axis equal, colormap(jet)
   %a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)

   hf=figure(10*Mac.plot_JM+3);
   pcolor(R(1:N,:)*Mac.R0EXP,Z(1:N,:)*Mac.R0EXP,real(J2(1:N,:))*fac_J), colorbar, shading interp,
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
   title('J^2 [MA/m^2]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   axis equal, colormap(jet)
   %a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)

   hf=figure(10*Mac.plot_JM+4);
   pcolor(R(1:N,:)*Mac.R0EXP,Z(1:N,:)*Mac.R0EXP,imag(J2(1:N,:))*fac_J), colorbar, shading interp,
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
   title('Im(J^2) [MA/m^2]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   axis equal, colormap(jet)
   %a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)

   hf=figure(10*Mac.plot_JM+5);
   pcolor(R(1:N,:)*Mac.R0EXP,Z(1:N,:)*Mac.R0EXP,real(J3(1:N,:))*fac_J), colorbar, shading interp,
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
   title('J^3 [MA/m^2]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   axis equal, colormap(jet)
   %a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)

   hf=figure(10*Mac.plot_JM+6);
   pcolor(R(1:N,:)*Mac.R0EXP,Z(1:N,:)*Mac.R0EXP,imag(J3(1:N,:))*fac_J), colorbar, shading interp,
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
   title('Im(J^3) [MA/m^2]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   axis equal, colormap(jet)
   %a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)

   hf=figure(10*Mac.plot_JM+7);
   pcolor(R(1:N,:)*Mac.R0EXP,Z(1:N,:)*Mac.R0EXP,real(JR(1:N,:))*fac_J), colorbar, shading interp,
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
   title('J_{RE} [MA/m^2]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   axis equal, colormap(jet)
   %a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)

   hf=figure(10*Mac.plot_JM+8);
   pcolor(R(1:N,:)*Mac.R0EXP,Z(1:N,:)*Mac.R0EXP,imag(JR(1:N,:))*fac_J), colorbar, shading interp,
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
   title('Im(J_{RE}) [MA/m^2]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   axis equal, colormap(jet)
   %a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)
 
   %plot surface averaged JRE
   if Mac.plot_Jpara
   hf=figure(10*Mac.plot_Jpara+7);
   %[a,II] = max(JRS);
   %JRS(1:II-1) = a + (a-JRS(1:II-1));
   plot(Mac.s(1:N),JRS*fac_J,[Mac.SS '-'],'LineWidth',2), hold on
   xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
   ylabel('<J_{||}> [MA/m^2]','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
   res = [Mac.s(1:N) JRS];
   save JparaSurfRE.txt res -ascii
   end
   end

   if 1==0
%  plot J1,J2,J3
   II=[1:Mac.Ns];
   MN = Mac.Mm; %MN = 2;  
   MM = MN - Mac.Mm(1) + 1;

   hf=figure(10*Mac.plot_JM+7);
   Y = real(JM1(II,MM)); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot((Mac.s(II(1:end-1))+Mac.s(II(2:end)))/2,Y(1:end-1,:),'LineWidth',1); hold on,
   ylabel('Re[J^1_{(m)}]','FontSize',18,'FontWeight','Bold')
   xlabel('s\equiv\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c) 
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end

   hf=figure(10*Mac.plot_JM+8);
   Y = imag(JM2(II,MM)); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(II),Y,'LineWidth',1); hold on,
   ylabel('Im[J^2_{(m)}]','FontSize',18,'FontWeight','Bold')
   xlabel('s\equiv\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c) 
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end

%  save m=2 harmonic of J2 
   %res = [Mac.s(II) Mac.q(II) real(JM2(II,MM)) imag(JM2(II,MM))];
   %save JM2.txt res -ascii


   hf=figure(10*Mac.plot_JM+9);
   Y = abs(JM3(II,MM)); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(II),Y,'LineWidth',1); hold on,
   ylabel('|J^3_{(m)}|','FontSize',18,'FontWeight','Bold')
   xlabel('s\equiv\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       %text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c) 
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   %axis([0.99 1.001 0 a(4)])

   end
end

if Mac.save_FEEDJ 
   [smin II] = min(abs(Mac.s-Mac.rw(1)));
   fid = fopen('FEEDJ','w');
   fprintf(fid,'%3i\n',Mac.Nm1);
   FEEDJ = [Mac.Mm'; real(JM2(II,:)); imag(JM2(II,:))];
   fprintf(fid,'%3i %16.8e %16.8e\n',FEEDJ);
   fclose(fid);

   data = [Mac.Mm real(JM2(II,:))' imag(JM2(II,:))' real(JM3(II,:))' imag(JM3(II,:))'];
   save JRW data -ascii 
end 

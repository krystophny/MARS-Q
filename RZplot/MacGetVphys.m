function [Vr,Vz,Vphi] = MacGetVphys(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,V1,V2,V3)

global Mac

N = Mac.Ns1;
R0 = R(1:N,:);  Z0 = Z(1:N,:);

% compute Vr,Vz,Vphi
jacobian(1,:) = jacobian(2,:);
%fac = Mac.B0EXP^2/Mac.R0EXP;
fac = 1.0;
Vr = (V1.*dRds(1:N,:) + V2.*dRdchi(1:N,:))*fac; Vr(1,:) = Vr(2,:);
Vz = (V1.*dZds(1:N,:) + V2.*dZdchi(1:N,:))*fac; Vz(1,:) = Vz(2,:);
Vphi = V3.*R0*fac;                     Vphi(1,:) = Vphi(2,:);

% plot real(Vr,Vz) at (R,Z) plane
if Mac.plot_V>0
   hf=figure(10*Mac.plot_V+0);
   pcolor(R0*Mac.R0EXP,Z0*Mac.R0EXP,real(Vz)), hold on
   shading interp, colormap(jet), colorbar
   axis equal
   xlabel('R [m]','FontSize',18,'FontWeight','Bold')
   ylabel('Z [m]','FontSize',18,'FontWeight','Bold')
   title('{\xi}_Z [mm]','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

   hf=figure(10*Mac.plot_V+1);
   II = find(Mac.s<=Mac.edge);  I=II(end); 
   plot(Mac.chi*180/pi,abs(Vz(I,:)),Mac.SS,'LineWidth',3), hold on
   xlabel('poloidal angle [deg]','FontSize',18,'FontWeight','Bold')
   ylabel('|{\xi}_Z| [mm]','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

   Rexpt = 1.94; %[m]
   II = find(Z0(end,:)>0); 
   [Y,I]=min(abs(R0(end,II)*Mac.R0EXP-Rexpt));
   Cexpt = Mac.chi(II(I))*180/pi;
   a=axis; plot([Cexpt Cexpt],[a(3) a(4)],'k--')

   figure(10*Mac.plot_V+0);
   R2 = R0(end,II(I))*Mac.R0EXP;
   Z2 = Z0(end,II(I))*Mac.R0EXP;
   plot(R2,Z2,'b+','MarkerSize',12)



   hf=figure(10*Mac.plot_V+2);
   pcolor(R0*Mac.R0EXP,Z0*Mac.R0EXP,real(Vr)), hold on
   shading interp, colormap(jet), colorbar
   axis equal
   xlabel('R [m]','FontSize',18,'FontWeight','Bold')
   ylabel('Z [m]','FontSize',18,'FontWeight','Bold')
   title('{\xi}_R [mm]','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

   hf=figure(10*Mac.plot_V+3);
   II = find(Mac.s<=Mac.edge);  I=II(end); 
   plot(Mac.chi*180/pi,real(Vr(I,:)),'b-','LineWidth',3), hold on
   plot(Mac.chi*180/pi,imag(Vr(I,:)),'b--','LineWidth',3), hold on
   xlabel('poloidal angle [deg]','FontSize',18,'FontWeight','Bold')
   ylabel('|{\xi}_R| [mm]','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

   Zexpt = 0; %[m]
   II = find(R0(end,:)>1); 
   [Y,I]=min(abs(Z0(end,II)*Mac.R0EXP-Zexpt));
   Cexpt = Mac.chi(II(I))*180/pi;
   a=axis; plot([Cexpt Cexpt],[a(3) a(4)],'k--')

   figure(10*Mac.plot_V+2);
   R2 = R0(end,II(I))*Mac.R0EXP;
   Z2 = Z0(end,II(I))*Mac.R0EXP;
   plot(R2,Z2,'b+','MarkerSize',12)

   hf=figure(10*Mac.plot_V + 7);
   % make original grid points and dense grid points to check differences
   phys_R = R0*Mac.R0EXP;
   phys_Z = Z0*Mac.R0EXP;
   phys_Vz = abs(Vz);
   cut_Z = linspace(0.65,0.8,401);
   cut_R = 1.94*ones(1,401);  % setting R=1.94m as Thomson axis
   cut_Vz = griddata(phys_R, phys_Z, phys_Vz, cut_R, cut_Z);

   plot(cut_Z, abs(cut_Vz), 'b-', 'LineWidth', 3), hold on
   title('{\xi}_Z [mm] vs Z at R=1.94m','FontSize',18,'FontWeight','Bold')
   xlabel('Z [m]','FontSize',18,'FontWeight','Bold')
   ylabel('|{\xi}_Z| [mm]','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

%  plot change of flux surfaces
   hf=figure(10*Mac.plot_V+8);
   sss = linspace(0,1,11); sss=sss(2:end);
   ccc = linspace(0,1,length(sss));
   for k=1:length(sss)
       [Y,II]=min(abs(Mac.s-sss(k)));
       col = [1-ccc(k) 0 ccc(k)];
       plot(R0(II,:)*Mac.R0EXP,Z0(II,:)*Mac.R0EXP,'--','Color',col,'LineWidth',2), hold on
       plot(R0(II,:)*Mac.R0EXP+real(Vr(II,:))/1000,Z0(II,:)*Mac.R0EXP+real(Vz(II,:))/1000,'-','Color',col,'LineWidth',2), hold on
   end     
   xlabel('R [m]','FontSize',18,'FontWeight','Bold')
   ylabel('Z [m]','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
   axis equal
end




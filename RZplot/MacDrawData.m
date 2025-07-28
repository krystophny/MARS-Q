% plot RMP comparison data (AUG vs. D3D)

dtyp = 3;

if dtyp == 1
%data='ProfQ';
data='ProfRot';
%data='B1m2D';
%data='B1m1D';
%data='Xn2D';
%data='X1m1D';
elseif dtyp == 2  
   data='B1res';
elseif dtyp == 3 
   data='XnSurfGeom';  
elseif dtyp == 4
   %data='B1m1D';
   %data='X1m1D';
   %data='ProfRot';
   %data='Xn2D';
    data='B1m2D';
elseif dtyp == 5  
   data='B1res';
elseif dtyp == 6 
   data='XnSurfGeom';  
end

cass='PLS_EF';
%cass='RMP_WF_Q95';

if dtyp == 1
   openfig(['../AUG30835/Data/' data '_AUG30835.3200_' cass '.fig'])
   openfig(['../AUG33133/Data/' data '_AUG33133.2600_' cass '.fig'])
   openfig(['../AUG33133/Data/' data '_AUG33133.3000_' cass '.fig'])
   openfig(['../D3D164362/Data/' data '_D3D164362.3400_' cass '.fig'])
   openfig(['../D3D164277/Data/' data '_D3D164277.2500_' cass '.fig'])
elseif dtyp == 2
   d1=load(['../AUG30835/Data/' data '_AUG30835.3200_' cass]);
   d2=load(['../AUG33133/Data/' data '_AUG33133.2600_' cass]);
   d3=load(['../AUG33133/Data_New/' data '_AUG33133.3000_' cass '_N']);
   d4=load(['../D3D164362/Data/' data '_D3D164362.3400_' cass]);
   d5=load(['../D3D164277/Data/' data '_D3D164277.2500_' cass]);

   hf=figure(1);
   plot(d1(:,2).^2,d1(:,6)*1e+4,'r-o','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','r'), hold on,
   %plot(d2(:,2).^2,d2(:,6)*1e+4,'b--s','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','b'), hold on,
   plot(d3(:,2).^2,d3(:,6)*1e+4,'k-.o','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','k'), hold on,
   %plot(d4(:,2).^2,d4(:,6)*1e+4,'k-d','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','k'), hold on,
   plot(d5(:,2).^2,d5(:,6)*1e+4,'b--s','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','b'), hold on,
   xlabel('\psi_p','FontSize',18,'FontWeight','Bold'),
   ylabel('b^1_{res} x 10^4','FontSize',18,'FontWeight','Bold'),
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),
   %hl=legend('AUG30835.3200','AUG33133.2600','AUG33133.3000','D3D164362.3400','D3D164277.2500');
   hl=legend('AUG30835.3200','AUG33133.3000','D3D164277.2500');
   set(hl,'FontSize',14),
elseif dtyp == 3 & cass(1:3)~='VAC'
   d1=load(['../AUG30835/Data/' data '_AUG30835.3200_' cass]);
   %d2=load(['../AUG33133/Data/' data '_AUG33133.2600_' cass]);
   d3=load(['../AUG33133/Data_New/' data '_AUG33133.3000_' cass '_N']);
   %d4=load(['../D3D164362/Data/' data '_D3D164362.3400_' cass]);
   %d5=load(['../D3D164277/Data/' data '_D3D164277.2500_' cass]);

   hf=figure(2);
   plot(d1(:,1),d1(:,2),'r-','LineWidth',4), hold on,
   %plot(d2(:,1),d2(:,2),'b--','LineWidth',3), hold on,
   plot(d3(:,1),d3(:,2),'k-.','LineWidth',4), hold on,
   %plot(d4(:,1),d4(:,2),'k-','LineWidth',3), hold on,
   %plot(d5(:,1),d5(:,2),'b--','LineWidth',3), hold on,
   xlabel('Poloidal angle [deg]','FontSize',18,'FontWeight','Bold'),
   ylabel('\xi_n [mm]','FontSize',18,'FontWeight','Bold'),
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),
   %hl=legend('AUG30835.3200','AUG33133.2600','AUG33133.3000','D3D164362.3400','D3D164277.2500');
   hl=legend('AUG30835.3200','AUG33133.3000');
   %set(hl,'FontSize',14),
elseif dtyp == 4
   aa=linspace(0,1,11);
   for k=11%1:length(aa)
       a = num2str(aa(k));   
       openfig(['../AUG30835/Data/' data '_AUG30835.3200_a=' a '_' cass '.fig'])
   end
elseif dtyp == 5
   hf=figure(1);
   aa=linspace(0,1,11);
   for k=1:length(aa)
       a = num2str(aa(k));   
       d=load(['../AUG30835/Data/' data '_AUG30835.3200_a=' a '_' cass]);
       plot(d(:,2).^2,d(:,6)*1e+4,'-s','LineWidth',2,'Color',[(11-k)/10 0 (k-1)/10],'MarkerSize',9,'MarkerFaceColor',[(11-k)/10 0 (k-1)/10]), hold on,
   end
   xlabel('\psi_p','FontSize',18,'FontWeight','Bold'),
   ylabel('b^1_{res} x 10^4','FontSize',18,'FontWeight','Bold'),
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),
   hl=legend('\alpha=0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1');
   set(hl,'FontSize',14),
elseif dtyp == 6
   aa=linspace(0,1,11);
   hf=figure(2);
   for k=1:length(aa)
       a = num2str(aa(k));   
       d=load(['../AUG30835/Data/' data '_AUG30835.3200_a=' a '_' cass]);
       plot(d(:,1),d(:,2),'-','LineWidth',2,'Color',[(11-k)/10 0 (k-1)/10]), hold on,
   end
   xlabel('geometric angle [deg]','FontSize',18,'FontWeight','Bold'),
   ylabel('\xi_n [mm]','FontSize',18,'FontWeight','Bold'),
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),
   hl=legend('\alpha=0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1');
   set(hl,'FontSize',14),
end

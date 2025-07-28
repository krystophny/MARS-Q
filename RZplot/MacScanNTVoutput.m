% post-process TORQUE_?? data: the net JxB,NTV and REY torques in [Nm]
% with 8x8 grid scan of coil phasing for ITER

IRMP   = 10; %[kAt]
kprint = 1;

SS = 'ABC';
%for k1=1:length(SS)
%for k2=1:3
%for k3=1:4
    
for k1=1:1
for k2=3:3
for k3=2:2
    
close all,

SDIR = ['/home/liuy/Work/ITER_NTV/' SS(k1) '/rot' int2str(k2) '/n' int2str(k3) '/'];

eval(['!sed -i s/''COMPUTE TOTAL JXB TORQUE:''/''' '''/g ' SDIR 'TORQUE_*']);
eval(['!sed -i s/''TOTAL NTV TORQUE:''/''' '''/g ' SDIR 'TORQUE_*']);
eval(['!sed -i s/''THERMAL ION NTV TORQUE:''/''' '''/g ' SDIR 'TORQUE_*']);
eval(['!sed -i s/''THERMAL ELECTRON NTV TORQUE:''/''' '''/g ' SDIR 'TORQUE_*']);
eval(['!sed -i s/''COMPUTE TOTAL REYNOLDS TORQUE:''/''' '''/g ' SDIR 'TORQUE_*']);

Tjxb  = zeros(9,9);
Trey  = zeros(9,9);
Tntvi = zeros(9,9);
Tntve = zeros(9,9);

phiU = linspace(0,360,9);
phiL = linspace(0,360,9);

for k=1:8
for l=1:8
    a = load([SDIR 'TORQUE_' int2str(k) int2str(l)]);
    Tjxb(k,l)  = a(1,2);
    Trey(k,l)  = a(6,2);
    Tntvi(k,l) = a(4,2);
    Tntve(k,l) = a(5,2);
end
end

Tjxb(:,9)  =  Tjxb(:,1); Tjxb(9,:)  =  Tjxb(1,:);
Trey(:,9)  =  Trey(:,1); Trey(9,:)  =  Trey(1,:);
Tntvi(:,9) =  Tntvi(:,1);Tntvi(9,:) =  Tntvi(1,:);
Tntve(:,9) =  Tntve(:,1);Tntve(9,:) =  Tntve(1,:);

Tjxb  = transpose(Tjxb)*IRMP^2;
Trey  = transpose(Trey)*IRMP^2;
Tntvi = transpose(Tntvi)*IRMP^2;
Tntve = transpose(Tntve)*IRMP^2;

hf = figure(1);
pcolor(phiU,phiL,Tjxb), 
colorbar, hold on, shading interp
axis([0 360 0 360]), 
xlabel('\Phi_U [deg.]','FontSize',18,'FontWeight','Bold')
ylabel('\Phi_L [deg.]','FontSize',18,'FontWeight','Bold')
title('JXB torque [Nm]','FontSize',14)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
if kprint>0, print(1,[SDIR 'torque_jxb'],'-djpeg'), end

hf = figure(2);
pcolor(phiU,phiL,Trey), 
colorbar, hold on, shading interp
axis([0 360 0 360]), 
xlabel('\Phi_U [deg.]','FontSize',18,'FontWeight','Bold')
ylabel('\Phi_L [deg.]','FontSize',18,'FontWeight','Bold')
title('REY torque [Nm]','FontSize',14)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
if kprint>0, print(2,[SDIR 'torque_rey'],'-djpeg'), end

hf = figure(3);
pcolor(phiU,phiL,Tntvi), 
colorbar, hold on, shading interp
axis([0 360 0 360]), 
xlabel('\Phi_U [deg.]','FontSize',18,'FontWeight','Bold')
ylabel('\Phi_L [deg.]','FontSize',18,'FontWeight','Bold')
title('i-NTV torque [Nm]','FontSize',14)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
if kprint>0, print(3,[SDIR 'torque_ntvi'],'-djpeg'), end

hf = figure(4);
pcolor(phiU,phiL,Tntve), 
colorbar, hold on, shading interp
axis([0 360 0 360]), 
xlabel('\Phi_U [deg.]','FontSize',18,'FontWeight','Bold')
ylabel('\Phi_L [deg.]','FontSize',18,'FontWeight','Bold')
title('e-NTV torque [Nm]','FontSize',14)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
if kprint>0, print(4,[SDIR 'torque_ntve'],'-djpeg'), end

end
end
end

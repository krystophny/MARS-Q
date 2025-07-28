function [Rc,Zc] = MacGetRZcoil(RM,ZM)

global Mac Acad

if Mac.Nm0<Mac.Nm2, Mac.Nm2 = Mac.Nm0; end

% Define the coils input for MARS-F
IFEED = [];
for k=1:size(Mac.coil,1)
  [smin,II] = min(abs(Mac.s-Mac.coil(k,1)));
  IFEED = [IFEED II-Mac.Ns1+1];
end
FCCHI = (Mac.coil(:,2) + Mac.coil(:,3))/2;
FWCHI = Mac.coil(:,3) - Mac.coil(:,2);

if Mac.resetCoil > 0
disp(['FCCHI=',num2str(FCCHI')])
disp(['FWCHI=',num2str(FWCHI')])
disp(['IFEED=',int2str(IFEED)])
disp(['SUMMARY=' num2str([FCCHI' FWCHI']) '  ' int2str(IFEED)])

Acad.FCCHI = FCCHI;
Acad.FWCHI = FWCHI;
Acad.IFEED = IFEED;
end

m = [0:1:Mac.Nm2-1]';
Nchi = 2;
Rc = zeros(size(Mac.coil,1),Nchi);
Zc = Rc;
for k=1:size(Mac.coil,1)
  rc   = Mac.coil(k,1);
  chi1 = Mac.coil(k,2)*pi;
  chi2 = Mac.coil(k,3)*pi;

  chi = linspace(chi1,chi2,Nchi);
  expmchi = exp(m*chi*i);

  [smin,II] = min(abs(Mac.s-rc));
  
  Rc(k,:) = real(RM(II,1:Mac.Nm2)*expmchi);
  Zc(k,:) = real(ZM(II,1:Mac.Nm2)*expmchi);
end

% plot RZ coordinates for coils
if Mac.plot_coil>0
   figure(Mac.plot_coil)

   if Mac.resetCoil
     RZ = Mac.coilN;
     for k=1:size(RZ,1);
          plot([RZ(k,1) RZ(k,3)],[RZ(k,2) RZ(k,4)],'r-s','LineWidth',1,'MarkerSize',9,'MarkerFaceColor','r'), hold on
     end
   end

   plot(Rc(:,1)*Mac.R0EXP,Zc(:,1)*Mac.R0EXP,'b+','LineWidth',2,'MarkerSize',9,'MarkerFaceColor','b'), hold on
   plot(Rc(:,end)*Mac.R0EXP,Zc(:,end)*Mac.R0EXP,'b+','LineWidth',2,'MarkerSize',9,'MarkerFaceColor','b'), hold on
   axis equal
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
end

if Mac.resetCoil > 0
CoilRZ = [Rc(:,1)*Mac.R0EXP Zc(:,1)*Mac.R0EXP
          Rc(:,end)*Mac.R0EXP Zc(:,end)*Mac.R0EXP]
end

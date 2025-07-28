function [Brho,Bchi,Bphi] = MacGetBphysT(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,B1,B2,B3)
%% Get B-field components in Toroidal coordinates

global Mac

sqrtG11 = sqrt(dRds.^2 + dZds.^2);
sqrtG22 = sqrt(dRdchi.^2 + dZdchi.^2);
sqrtG33 = R;
Brho = B1./sqrtG22./R;  Brho(1,:) = Brho(2,:);
Bchi = B2./sqrtG11./R;  Bchi(1,:) = Bchi(2,:);
Bphi = B3.*sqrtG33./jacobian;  Bphi(1,:) = Bphi(2,:);

% plot real part of Brho
if Mac.plot_B > 0
   figure(Mac.plot_B)
   subplot(1,2,1)
   B = real(Brho);
   Bmax = max(max(abs(Brho)));
   B = B/Bmax;
   pcolor(R*Mac.R0EXP,Z*Mac.R0EXP,B), hold on,
   axis equal
   axis([0 5.0 -3.3 3.3])
   colorbar
   shading interp
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
   title('Re[B_r]','FontSize',14)

   subplot(1,2,2)
   B = imag(Brho);
   Bmax = max(max(abs(Brho)));
   B = B/Bmax;
   pcolor(R*Mac.R0EXP,Z*Mac.R0EXP,B), hold on,
   axis equal
   axis([0 5.0 -3.3 3.3])
   colorbar
   shading interp
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
   title('Im[B_r]','FontSize',14)
end




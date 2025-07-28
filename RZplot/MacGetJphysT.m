function [Jrho,Jchi,Jphi] = MacGetJphysT(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,J1,J2,J3)
%% Get J-field components in Toroidal coordinates

global Mac

sqrtG11 = sqrt(dRds.^2 + dZds.^2);
sqrtG22 = sqrt(dRdchi.^2 + dZdchi.^2);
sqrtG33 = R;

N = size(J1,1);
Jrho = J1.*sqrtG11(1:N,:)./jacobian(1:N,:);  Jrho(1,:) = Jrho(2,:);
Jchi = J2.*sqrtG22(1:N,:)./jacobian(1:N,:);  Jchi(1,:) = Jchi(2,:);
Jphi = J3.*sqrtG33(1:N,:)./jacobian(1:N,:);  Jphi(1,:) = Jphi(2,:);

% plot real part of J at (chi,phi) plane
if Mac.plot_JW2 > 0
   figure(Mac.plot_JW2)
   nrw = length(Mac.rw);
   wt = 0;
   [pp,cc] = meshgrid(Mac.phi,Mac.chi);
   pp = pp*180/pi;
   cc = cc*180/pi;
   for k=1:nrw
   subplot(nrw,1,k)
   [smin,II] = min(abs(Mac.s-Mac.rw(k)));
   Jwc = abs(exp(i*Mac.n*Mac.phi+i*wt)*Jchi(II,:))';
   Jwp = abs(exp(i*Mac.n*Mac.phi+i*wt)*Jphi(II,:))';
   Jwa = sqrt(Jwp.^2+Jwc.^2);
   Jmax= max(max(Jwa));

   %quiver(pp,cc,Jwp./Jwa,Jwc./Jwa,0.8), hold on
   pcolor(pp,cc,Jwa/Jmax), shading interp, hold on,
   colorbar,  colormap(hot)

   xlabel('geometric toroidal angle','FontSize',16,'FontWeight','Bold')
   ylabel('geometric poloidal angle','FontSize',16,'FontWeight','Bold')
   %title('inner wall','FontSize',14)

   %find poloidal angles for top,bottom,LFS,HFS
   [Y,k]=min(R(II,:)); chi_HFS=Mac.chi(k)*180/pi;
   [Y,k]=max(R(II,:)); chi_LFS=Mac.chi(k)*180/pi;
   [Y,k]=min(Z(II,:)); chi_bot=Mac.chi(k)*180/pi;
   [Y,k]=max(Z(II,:)); chi_top=Mac.chi(k)*180/pi;
   axis([0 400 -180 180])
   text(365,chi_HFS,'HFS','FontSize',16,'FontWeight','Bold')
   text(365,chi_LFS,'LFS','FontSize',16,'FontWeight','Bold')
   text(365,chi_bot,'BOT','FontSize',16,'FontWeight','Bold')
   text(365,chi_top,'TOP','FontSize',16,'FontWeight','Bold')
   end
end

% plot Jphi along the outboard mid-plane
if Mac.plot_Jphi > 0
   hf=figure(10*Mac.plot_Jphi+0);
   [X,I2]= min(abs(Mac.chi));
   I1    = [1:Mac.Ns1-1];
   %plot(Mac.s(I1),real(Jphi(I1,I2)),'b-','LineWidth',2), hold on,
   %plot(Mac.s(I1),imag(Jphi(I1,I2)),'k-','LineWidth',2), hold on,
   plot(Mac.s(I1), abs(Jphi(I1,I2)),'r-','LineWidth',2), hold on,
   xlabel('s','FontSize',16)
   ylabel('J_{\phi}','FontSize',16)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
end




  data = load([SDIR 'XnSurf.txt']);
  chi  = data(:,1); 
  Rn   = data(:,2);
  Zn   = data(:,3);
  Xn   = data(:,4)+data(:,5)*i;

  %normal displacement at the plasma surface, along poloidal angle
  hf=figure(47);
  VV = abs(Xn);
  x = chi*180/pi;
  plot(x,VV*1000,'-','LineWidth',2,'Color',c), hold on,
  xlabel('equal-arc. poloidal angle [deg.]','FontSize',16,'FontWeight','Bold')
  ylabel('|\xi_n| [mm/kA]','FontSize',16,'FontWeight','Bold')
  ha = get(hf,'CurrentAxes');
  set(ha,'FontSize',14,'FontWeight','Bold')

  %normal displacement along the plasma surface, at R-Z plane
  hf=figure(48);
  fac_scale = 100;
  VV = imag((Xn(1:end-1)+Xn(2:end)))/2;
  Y1 = Rn(1:end-1) + Zn(1:end-1)*i;
  Y2 = Rn(2:end)   + Zn(2:end)*i;
  Y0 = (Y1+Y2)/2;
  YD = Y2-Y1;
  YY = Y0 + VV.*YD./abs(YD)*(-i)*fac_scale;
  plot(real(Y0),imag(Y0),'k--','LineWidth',0.5), hold on,
  plot(real(YY),imag(YY),'-','LineWidth',2,'Color',c), hold on,
  xlabel('R [m]','FontSize',16,'FontWeight','Bold')
  ylabel('Z [m]','FontSize',16,'FontWeight','Bold')
  ha = get(hf,'CurrentAxes');
  set(ha,'FontSize',14,'FontWeight','Bold')
  text(1.25,-1,'\xi_n x100','FontSize',18,'FontWeight','Bold')
  axis equal
  axis([0 2 -1.25 1.25]) 

function [Bmid,Btmid,BBtave,Btave] = MacGetNTVBeq(R,Z,dRdchi,dZdchi,jacobian)

global Mac

G22 = dRdchi(1:Mac.Ns1,:).^2 + dZdchi(1:Mac.Ns1,:).^2;
Jac = jacobian(1:Mac.Ns1,:);
RR  = R(1:Mac.Ns1,:);

%compute equilibrium fields magnitude on (R,Z) 
onevec = ones(1,size(RR,2));
BB = sqrt(Mac.DPSIDS.^2*onevec.*G22./Jac.^2 + Mac.T.^2*onevec./RR.^2);
BB(1,:) = sqrt(Mac.DPSIDS(2).^2*onevec.*G22(1,:)./Jac(2,:).^2 + Mac.T(1).^2*onevec./RR(1,:).^2);
BB(1,:) = sum(BB(1,:))/size(BB,2);
Bt = Mac.T*onevec./RR;

%compute equilibrium fields magnitude along R at the midplane Z=Zc
Rc = R(1,1); Zc = Z(1,1);
Bmid  = zeros(Mac.Ns1,1); Bmid(1)  = sum(BB(1,:))/size(BB,2);
Btmid = zeros(Mac.Ns1,1); Btmid(1) = Bt(1,1);

for k=2:Mac.Ns1
  JJ = find(R(k,:) >= Rc);
  yy = BB(k,JJ);
  yt = Bt(k,JJ);
  ZZ = Z(k,JJ);
  
  JJ = find( (ZZ(1:end-1)-Zc).*(ZZ(2:end)-Zc) <= 0 );
  Bmid(k)  = yy(JJ) + (Zc-ZZ(JJ))*(yy(JJ+1)-yy(JJ))/(ZZ(JJ+1)-ZZ(JJ));
  Btmid(k) = yt(JJ) + (Zc-ZZ(JJ))*(yt(JJ+1)-yt(JJ))/(ZZ(JJ+1)-ZZ(JJ));
end

%compute surface-averaged equilibrium fields
%notice the weighting kernal J/R due to coordinates transform
hh = Mac.chi(2) - Mac.chi(1);
yy = jacobian(2:Mac.Ns1,:)./BB(2:Mac.Ns1,:)./Bt(2:Mac.Ns1,:)./R(2:Mac.Ns1,:);
yt = jacobian(2:Mac.Ns1,:)./Bt(2:Mac.Ns1,:)./R(2:Mac.Ns1,:);
yd = jacobian(2:Mac.Ns1,:)./R(2:Mac.Ns1,:);
zy = hh*(sum(yy,2)-yy(:,1)*0.5-yy(:,end)*0.5);
zt = hh*(sum(yt,2)-yt(:,1)*0.5-yt(:,end)*0.5);
zd = hh*(sum(yd,2)-yd(:,1)*0.5-yd(:,end)*0.5);
BBtave = [1.0/BB(1,1)/Bt(1,1); zy./zd];
Btave  = [1.0/Bt(1,1); zt./zd];

%re-normalization
Bmid   = Bmid*Mac.B0EXP;
Btmid  = Btmid*Mac.B0EXP;
BBtave = BBtave/Mac.B0EXP^2;
Btave  = Btave/Mac.B0EXP;

if Mac.plot_NTVBeq > 0
  figure(Mac.plot_NTVBeq)
  plot(Mac.s(1:Mac.Ns1),Bmid,'b-','LineWidth',2), hold on,
  plot(Mac.s(1:Mac.Ns1),BBtave,'r-','LineWidth',2), hold on,
  plot(Mac.s(1:Mac.Ns1),Btmid,'b--','LineWidth',2), hold on,
  plot(Mac.s(1:Mac.Ns1),Btave,'r--','LineWidth',2), hold on,
  xlabel('s','FontSize',14)
  ylabel('Bmid(blue) & Bave(red)','FontSize',14)
end



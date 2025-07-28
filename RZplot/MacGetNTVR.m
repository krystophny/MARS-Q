function [Rmid,Rave] = MacGetNTVR(R,Z,jacobian)

global Mac

%compute R at mid-plane as function of Mac.s
Rc = R(1,1); Zc = Z(1,1);
Rmid = zeros(Mac.Ns1,1);
Rmid(1) = Rc;

for k=2:Mac.Ns1
  JJ = find(R(k,:) >= Rc);
  RR = R(k,JJ);
  ZZ = Z(k,JJ);
  
  JJ = find( (ZZ(1:end-1)-Zc).*(ZZ(2:end)-Zc) <= 0 );
  R0 = RR(JJ) + (Zc-ZZ(JJ))*(RR(JJ+1)-RR(JJ))/(ZZ(JJ+1)-ZZ(JJ));
  Rmid(k) = R0;
end

%compute surface-averaged 1/R^2
%notice the weighting kernal J/R due to coordinates transform
hh = Mac.chi(2) - Mac.chi(1);
yy = jacobian(2:Mac.Ns1,:)./R(2:Mac.Ns1,:).^3; 
zn = hh*(sum(yy,2)-yy(:,1)*0.5-yy(:,end)*0.5); 
yy = jacobian(2:Mac.Ns1,:)./R(2:Mac.Ns1,:); 
zd = hh*(sum(yy,2)-yy(:,1)*0.5-yy(:,end)*0.5); 
Rave = [1/Rc^2; zn./zd];

%re-normalization
Rmid = Rmid*Mac.R0EXP;
Rave = Rave/Mac.R0EXP^2;


if Mac.plot_NTVR > 0
  figure(Mac.plot_NTVR)
  plot(Mac.s(1:Mac.Ns1),Rmid,'b-','LineWidth',2), hold on,
  plot(Mac.s(1:Mac.Ns1),Rave,'r-','LineWidth',2), hold on,
  xlabel('s','FontSize',14)
  ylabel('Rmid(blue) & Rave(red)','FontSize',14)
end

%data = [Mac.s(1:Mac.Ns1) Rmid];
%save dataSR data -ascii -double 

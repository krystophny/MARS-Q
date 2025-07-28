function [Jr,Jz,Jphi] = MacGetJphysC(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,J1,J2,J3)

global Mac

N = Mac.Ns;
R0 = R(1:N,:);  Z0 = Z(1:N,:);

% compute Jr,Jz,Jphi, jacobian from (R,Z,phi) --> (s,chi,phi)
jacobian(1,:) = jacobian(2,:);
Jr = (J1.*dRds(1:N,:) + J2.*dRdchi(1:N,:))./jacobian(1:N,:);  Jr(1,:) = Jr(2,:);
Jz = (J1.*dZds(1:N,:) + J2.*dZdchi(1:N,:))./jacobian(1:N,:);  Jz(1,:) = Jz(2,:);
Jphi = J3.*R0./jacobian(1:N,:);                     Jphi(1,:) = Jphi(2,:);

% plot real(Jr,Jz) at (R,Z) plane
if Mac.plot_JP>0
   h = 3.0/Mac.Ns1;
   Jrr = real(Jr); Jzr = real(Jz);
   Jt = sqrt(Jrr.^2 + Jzr.^2);  
   %Jt = max(max(Jt));
   [II,JJ]=find(Jt==0.0);  Jt(II,JJ) = 1.0;
   R1 = R0 + h*Jrr./Jt;
   Z1 = Z0 + h*Jzr./Jt;

   figure(Mac.plot_JP)
   R2 = [R0(:) R1(:)]';
   Z2 = [Z0(:) Z1(:)]';

   plot(R0*Mac.R0EXP,Z0*Mac.R0EXP,'c-'), hold on,
   plot(R0'*Mac.R0EXP,Z0'*Mac.R0EXP,'y-'), hold on,
   plot(R2*Mac.R0EXP,Z2*Mac.R0EXP,'b-'), hold on,
   axis equal
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
end




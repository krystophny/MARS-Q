function Torq = MacGetNTVTorq(BMH,Rmid,Rave,Vti,Bmid,Btmid,BBtave,Btave,Wnm,nui,rot,peqi,Sv);

global Mac

%NTV force in plateau and collisional regimes
mups1 = 1.365;
C = sqrt(pi).*peqi./Vti.*rot.*Mac.q.*Rmid.^2.*Rave.*Bmid.*Btmid.*BBtave;
D = 0;
for m=1:size(BMH,2)
  mm = Mac.Mm(m);
  D = D + Mac.n^2*abs(BMH(:,m)).^2*mups1./(2*sqrt(pi)/3*nui.*Rmid.*Mac.q./Vti + mups1*abs(mm+Mac.n*Mac.q));
end
Fp = C.*D.*Sv./Bmid.^2; 

%NTV force in collisionless regime (due to trapped particles)
lamb1i = 13.708;
C = Btmid.*Rmid.*Btave.*Rave*lamb1i.*peqi/pi^1.5./nui.*rot.*(Mac.s(1:Mac.Ns1)/Mac.R0EXP).^1.5/sqrt(2)*Mac.n^2;
D = abs(BMH).^2.*Wnm;
Fnu = C.*sum(D,2).*Sv./Bmid.^2; 

II = 50;
res = [Btmid(II); Rmid(II); Btave(II); Rave(II); lamb1i; peqi(II); pi^1.5; nui(II); rot(II); (Mac.s(II)/Mac.R0EXP)^1.5; sqrt(2); Mac.n; sum(abs(BMH(II,:)).^2); sum(D(II,:)); Sv(II)]
 
%total torque
Torq = Rmid.*(Fp + Fnu);

if Mac.plot_NTV > 0
  figure(Mac.plot_NTV)
  plot(Rmid,Torq,'b-','LineWidth',2), hold on,
  %plot(Rmid,Rmid.*Fp,'b--','LineWidth',1), hold on,
  %plot(Rmid,Rmid.*Fnu,'b-.','LineWidth',1), hold on,
  xlabel('R [m]','FontSize',16),
  ylabel('torque density [N]','FontSize',16),
end


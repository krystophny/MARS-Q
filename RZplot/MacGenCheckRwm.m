%mm=[-9:39];

MacMainITER_Srwm
Bm = [mm' real(Bnm(:)) imag(Bnm(:)) real(Btm(:)) imag(Btm(:)) real(Bpm(:)) imag(Bpm(:))];
save BmTotal.asc Bm -ascii -double
BnmT = Bnm;  BtmT = Btm;  BpmT = Bpm; 

MacMainITER_Srbio
Bm = [mm' real(Bnm(:)) imag(Bnm(:)) real(Btm(:)) imag(Btm(:)) real(Bpm(:)) imag(Bpm(:))];
save BmPlasma.asc Bm -ascii -double
BnmP = Bnm;  BtmP = Btm;  BpmP = Bpm; 

MacMainITER_Swall
Bm = [mm' real(Bnm(:)) imag(Bnm(:)) real(Btm(:)) imag(Btm(:)) real(Bpm(:)) imag(Bpm(:))];
save BmWall.asc Bm -ascii -double
BnmW = Bnm;  BtmW = Btm;  BpmW = Bpm; 

BnmT2 = BnmP + BnmW;
BtmT2 = BtmP + BtmW;
BpmT2 = BpmP + BpmW;

resn = [mm' BnmT(:) BnmT2(:) abs(BnmT(:)-BnmT2(:))./abs(BnmT(:))];
rest = [mm' BtmT(:) BtmT2(:) abs(BtmT(:)-BtmT2(:))./abs(BtmT(:))];
resp = [mm' BpmT(:) BpmT2(:) abs(BpmT(:)-BpmT2(:))./abs(BpmT(:))];

% get B vs. theta angle
t = linspace(-pi,pi,937)';
expmt = exp(i*t*mm);
BnT  = expmt*BnmT(:); BtT  = expmt*BtmT(:); BpT  = expmt*BpmT(:);
BnP  = expmt*BnmP(:); BtP  = expmt*BtmP(:); BpP  = expmt*BpmP(:);
BnW  = expmt*BnmW(:); BtW  = expmt*BtmW(:); BpW  = expmt*BpmW(:);

BnT2 = BnP + BnW;
BtT2 = BtP + BtW;
BpT2 = BpP + BpW;

RelError = [sqrt(sum(abs(BnT-BnT2).^2))/sqrt(sum(abs(BnT).^2)) ...
            sqrt(sum(abs(BtT-BtT2).^2))/sqrt(sum(abs(BtT).^2)) ...
            sqrt(sum(abs(BpT-BpT2).^2))/sqrt(sum(abs(BpT).^2))];
disp(['MacGenCheckRwm: RelError=' num2str(RelError)])

figure
subplot(3,2,1), plot(t,real(BnT),'r-','LineWidth',2), hold on,
                plot(t,real(BnP),'b-','LineWidth',1), hold on,
                plot(t,real(BnW),'k-','LineWidth',1), hold on,
                plot(t,real(BnT2),'r--','LineWidth',1), hold on,
		ylabel('Re(Bn)','FontSize',14),
subplot(3,2,2), plot(t,imag(BnT),'r-','LineWidth',2), hold on,
                plot(t,imag(BnP),'b-','LineWidth',1), hold on,
                plot(t,imag(BnW),'k-','LineWidth',1), hold on,
                plot(t,imag(BnT2),'r--','LineWidth',1), hold on,
		ylabel('Im(Bn)','FontSize',14),
subplot(3,2,3), plot(t,real(BtT),'r-','LineWidth',2), hold on,
                plot(t,real(BtP),'b-','LineWidth',1), hold on,
                plot(t,real(BtW),'k-','LineWidth',1), hold on,
                plot(t,real(BtT2),'r--','LineWidth',1), hold on,
		ylabel('Re(Bt)','FontSize',14),
subplot(3,2,4), plot(t,imag(BtT),'r-','LineWidth',2), hold on,
                plot(t,imag(BtP),'b-','LineWidth',1), hold on,
                plot(t,imag(BtW),'k-','LineWidth',1), hold on,
                plot(t,imag(BtT2),'r--','LineWidth',1), hold on,
		ylabel('Im(Bt)','FontSize',14),
subplot(3,2,5), plot(t,real(BpT),'r-','LineWidth',2), hold on,
                plot(t,real(BpP),'b-','LineWidth',1), hold on,
                plot(t,real(BpW),'k-','LineWidth',1), hold on,
                plot(t,real(BpT2),'r--','LineWidth',1), hold on,
		ylabel('Re(Bp)','FontSize',14),
		xlabel('\theta','FontSize',14),
subplot(3,2,6), plot(t,imag(BpT),'r-','LineWidth',2), hold on,
                plot(t,imag(BpP),'b-','LineWidth',1), hold on,
                plot(t,imag(BpW),'k-','LineWidth',1), hold on,
                plot(t,imag(BpT2),'r--','LineWidth',1), hold on,
		ylabel('Im(Bp)','FontSize',14),
		xlabel('\theta','FontSize',14),

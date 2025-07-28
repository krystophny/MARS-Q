%SDIR = '/home/elf/elfliu/FUSION/Machines/ITER/Result/Scen4_V8/';
SSDIR = '';

eval(['load ' SSDIR 'ResMatrixTot_t.asc']);
A = ResMatrixTot_t;
eval(['load ' SSDIR 'ResMatrixTot_p.asc']);
B = ResMatrixTot_p;

M = size(A,1);
A = 1*A(:,1:M) + i*A(:,M+1:end);  
B = 1*B(:,1:M) + i*B(:,M+1:end);

SDIR = [SDIR0 'BPLASMA_Srwm/'];
MacMainITER_Srwm

Btm2 = A*Bnm(:);
Bpm2 = B*Bnm(:);

%rest = [Btm(:) Btm2 abs(Btm(:)-Btm2)./abs(Btm(:))]
%resp = [Bpm(:) Bpm2 abs(Bpm(:)-Bpm2)./abs(Bpm(:))]


% plot Bt and Bt2 vs. theta angle
t = linspace(-pi,pi,257)';
expmt = exp(i*t*mm);
Bn  = expmt*Bnm(:);
Bt  = expmt*Btm(:);
Bp  = expmt*Bpm(:);
Bt2 = expmt*Btm2(:);
Bp2 = expmt*Bpm2(:);


RelError_t = sqrt(sum(abs(Bt-Bt2).^2))/sqrt(sum(abs(Bt).^2));
RelError_p = sqrt(sum(abs(Bp-Bp2).^2))/sqrt(sum(abs(Bp).^2));
disp(['MacGenCheckTot: RelError=' num2str(RelError_t) ', ' num2str(RelError_p)])

figure
subplot(3,2,1), plot(t,real(Bn),'r-','LineWidth',2), hold on,
		ylabel('Re(Bn)','FontSize',14),
subplot(3,2,2), plot(t,imag(Bn),'r-','LineWidth',2), hold on,
		ylabel('Im(Bn)','FontSize',14),
subplot(3,2,3), plot(t,real(Bt),'r-',t,real(Bt2),'b--','LineWidth',2), hold on,
		ylabel('Re(Bt)','FontSize',14),
subplot(3,2,4), plot(t,imag(Bt),'r-',t,imag(Bt2),'b--','LineWidth',2), hold on,
		ylabel('Im(Bt)','FontSize',14),
subplot(3,2,5), plot(t,real(Bp),'r-',t,real(Bp2),'b--','LineWidth',2), hold on,
		ylabel('Re(Bp)','FontSize',14),
		xlabel('\theta','FontSize',14),
subplot(3,2,6), plot(t,imag(Bp),'r-',t,imag(Bp2),'b--','LineWidth',2), hold on,
		ylabel('Im(Bp)','FontSize',14),
		xlabel('\theta','FontSize',14),

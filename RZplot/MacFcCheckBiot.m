%SDIR = '/home/elf/elfliu/FUSION/Machines/ITER/Result/Scen4_V8/';
SSDIR = '';

eval(['load ' SSDIR 'ResMatrixBiot_n.asc']);
A = ResMatrixBiot_n;

M = size(A,1);
A = 1*A(:,1:M) + i*A(:,M+1:end); 

SDIR = [SDIR0 'BPLASMA_Srwm/'];
MacMainITER_Srwm
Bnm2 = A*Bnm(:);

MacMainITER_Srbio
resn = [Bnm(:) Bnm2 abs(Bnm(:)-Bnm2)./abs(Bnm(:))];

% plot Bn and Bn2 vs. theta angle
t = linspace(-pi,pi,257)';
expmt = exp(i*t*mm);
Bn  = expmt*Bnm(:);
Bn2 = expmt*Bnm2(:);

RelError = sqrt(sum(abs(Bn-Bn2).^2))/sqrt(sum(abs(Bn).^2));
disp(['MacGenCheckBiot: RelError=' num2str(RelError)])

figure
subplot(2,1,1), plot(t,real(Bn),'r-',t,real(Bn2),'b--','LineWidth',2), hold on,
		ylabel('Re(Bn)','FontSize',14),
subplot(2,1,2), plot(t,imag(Bn),'r-',t,imag(Bn2),'b--','LineWidth',2), hold on,
		ylabel('Im(Bn)','FontSize',14),
		xlabel('\theta','FontSize',14),

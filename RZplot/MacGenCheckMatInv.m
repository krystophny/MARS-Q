%SDIR = '/home/elf/elfliu/FUSION/Machines/ITER/Result/Scen4_V8/';
SSDIR = '';

eval(['load ' SSDIR 'ResMatrixBiot_n.asc']);

qq=load([SSDIR 'BmTotal.asc']);
bntot_tilde=qq(:,2)+j*qq(:,3);

qq=load([SSDIR 'BmPlasma.asc']);
bnex_tilde=qq(:,2)+j*qq(:,3);

qq=load([SSDIR 'BmWall.asc']);
bned_tilde=qq(:,2)+j*qq(:,3);

qq=load([SSDIR 'ResMatrixBiot_n.asc']);
Bn_p_tilde_marsf=qq(:,1:MGG)+j*qq(:,MGG+1:end);
Bn_p_tilde_marsf=Bn_p_tilde_marsf*1.00;
dB_En_dB_Tn_marsf=eye(size(Bn_p_tilde_marsf))-Bn_p_tilde_marsf;

bnp_test_marsf=Bn_p_tilde_marsf*bntot_tilde;

bntot_test_marsf=inv(dB_En_dB_Tn_marsf)*bned_tilde;

%mm = [-9:39];

t = linspace(-pi,pi,257)';
expmt = exp(i*t*mm);

figure
ba = bnex_tilde;
bb = bnp_test_marsf;
bar= expmt*ba(:);
bbr= expmt*bb(:);
subplot(3,1,1), plot(mm,abs(ba),'b-'), hold on,
                plot(mm,abs(bb),'r--'), hold on,
subplot(3,1,2), plot(t,real(bar),'b-'), hold on,
                plot(t,real(bbr),'r--'), hold on,
subplot(3,1,3), plot(t,imag(bar),'b-'), hold on,
                plot(t,imag(bbr),'r--'), hold on,
RelError = sqrt(sum(abs(bar-bbr).^2))/sqrt(sum(abs(bar).^2));
disp(['MacGenCheckMatInv:bnex: RelError=' num2str(RelError)])

figure
ba = bntot_tilde;
bb = bntot_test_marsf;
bar= expmt*ba(:);
bbr= expmt*bb(:);
subplot(3,1,1), plot(mm,abs(ba),'b-'), hold on,
                plot(mm,abs(bb),'r--'), hold on,
subplot(3,1,2), plot(t,real(bar),'b-'), hold on,
                plot(t,real(bbr),'r--'), hold on,
subplot(3,1,3), plot(t,imag(bar),'b-'), hold on,
                plot(t,imag(bbr),'r--'), hold on,
RelError = sqrt(sum(abs(bar-bbr).^2))/sqrt(sum(abs(bar).^2));
disp(['MacGenCheckMatInv:bntot: RelError=' num2str(RelError)])

figure
ba = bned_tilde;
bb = bntot_tilde-bnp_test_marsf;
bar= expmt*ba(:);
bbr= expmt*bb(:);
subplot(3,1,1), plot(mm,abs(ba),'b-'), hold on,
                plot(mm,abs(bb),'r--'), hold on,
                plot(mm,abs(bntot_tilde-bnex_tilde),'b--'), hold on,
subplot(3,1,2), plot(t,real(bar),'b-'), hold on,
                plot(t,real(bbr),'r--'), hold on,
subplot(3,1,3), plot(t,imag(bar),'b-'), hold on,
                plot(t,imag(bbr),'r--'), hold on,
RelError = sqrt(sum(abs(bar-bbr).^2))/sqrt(sum(abs(bar).^2));
disp(['MacGenCheckMatInv:bned: RelError=' num2str(RelError)])


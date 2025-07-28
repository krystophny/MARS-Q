% generate and test plasma response matrices

global Mac

mm   = [-29:29];
tt   = linspace(-pi,pi,201)';
MGG  = length(mm);

IB1 = 264;  %always associated with data B123U_*
IB2 = 63;   %always associated with data B123U2_*

Fexp = exp(i*tt*mm);
ht   = tt(2)-tt(1);

MacMainD3D141069
Mac.Mm = mm';
[T,IIT]=MacPrmGetT(R,Z,IB1);

SDIR = '/.automount/funsrv1/root/home/yliu/D3D141069/Rfa/';
eval(['load ' SDIR 'B123U_ALLv'])
eval(['load ' SDIR 'B123U_ALLp'])
eval(['load ' SDIR 'B123U_Cv'])
eval(['load ' SDIR 'B123U_Cp'])
eval(['load ' SDIR 'B123U_Iv'])
eval(['load ' SDIR 'B123U_Ip'])

% generate response matrix MPRM
BGGv = zeros(MGG,MGG); BGGp = BGGv;
for GenNo=1:MGG
    BGGv(:,GenNo) = MacPrmBr(B123U_ALLv(GenNo,:),IB1,R,dRdchi,dZdchi,T,IIT);     
    BGGp(:,GenNo) = MacPrmBr(B123U_ALLp(GenNo,:),IB1,R,dRdchi,dZdchi,T,IIT);     
end

condnum_BGGv = cond(BGGv)

MPRM = BGGp * inv(BGGv);

% svd decomposition
[U,S,V] = svd(MPRM);

hf=figure(10);
semilogy(diag(S),'ro','LineWidth',2,'MarkerSize',7), hold on,
xlabel('mode number','FontSize',16,'FontWeight','Bold')
ylabel('SVD singular value','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

% compute the radial component of the vacuum n=1 error field
% and test the true error field response
BEvF = MacPrmVef(R,Z,IB1,T,IIT);
BEpnF = MPRM * BEvF;

BEvR = Fexp*BEvF;
BEpnR = Fexp*BEpnF;

hf=figure(20);
CEv = V' * BEvF;  
plot(abs(CEv),'r-o','LineWidth',2,'MarkerSize',7), hold on,
xlabel('mode number','FontSize',16,'FontWeight','Bold')
ylabel('mode amplitude','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf=figure(21);
plot(tt*180/pi,abs(BEvR),'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(BEpnR),'r-','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('|b_n(n=1)| [Gauss]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('vacuum field','total response field')

% test the true C-coil response
BCvF  = MacPrmBr(B123U_Cv,IB1,R,dRdchi,dZdchi,T,IIT);
BCpF  = MacPrmBr(B123U_Cp,IB1,R,dRdchi,dZdchi,T,IIT);

BCpnF = MPRM * BCvF;

BCvR  = Fexp*BCvF;
BCpR  = Fexp*BCpF;
BCpnR = Fexp*BCpnF;

hf=figure(30);
CCv = V' * BCvF;  
plot(abs(CCv),'r-o','LineWidth',2,'MarkerSize',7), hold on,
xlabel('mode number','FontSize',16,'FontWeight','Bold')
ylabel('mode amplitude','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf=figure(31);
plot(tt*180/pi,abs(BCvR),'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(BCpnR),'r-','LineWidth',2), hold on,
%plot(tt,abs(BCpR),'r--','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('|b_n(n=1)| [Gauss/kA-turn]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('vacuum field','total response field')

d = -CEv(1)/CCv(1);
SOL_C = log10(sqrt(sum(abs(d*BCpnR+BEpnR).^2)*ht));
d = [abs(d) angle(d)*180/pi d/abs(d) SOL_C];
SOLUTION_C = d 
hf=figure(32);
plot(tt*180/pi,abs(d(1)*BCvR*d(3)+BEvR),'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(d(1)*BCpnR*d(3)+BEpnR),'r-','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('|b_n(n=1)| [Gauss]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('vacuum field','total response field')

A       = linspace(0,4,101);
P       = linspace(0,360,37);
[AA,PP] = meshgrid(A,P);
PPD     = exp(i*PP*pi/180);
SOL     = AA;
for k=1:size(AA,1)
    for l=1:size(AA,2)
        SOL(k,l) = sum(abs(AA(k,l)*BCpnR*PPD(k,l)+BEpnR).^2);
    end
end
SOL = log10(sqrt(SOL*ht));

[Y,I] = min(SOL,[],2); [Y,K] = min(Y);
d = [AA(K,I(K)) PP(K,I(K)) PPD(K,I(K)) SOL(K,I(K))]; 
OPT_SOLUTION_C = d

hf=figure(33);
pcolor(AA,PP,SOL), hold on, shading interp
contour(AA,PP,SOL,20,'b-'), hold on,
xlabel('C-coil amplitude [kA-turn]','FontSize',16,'FontWeight','Bold')
ylabel('C-coil Phase [degree]','FontSize',16,'FontWeight','Bold')
colorbar,  colormap(hot)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf=figure(34);
plot(tt*180/pi,abs(d(1)*BCvR*d(3)+BEvR),'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(d(1)*BCpnR*d(3)+BEpnR),'r-','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('|b_n(n=1)| [Gauss]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('vacuum field','total response field')


% test the true I-coil response
BIvF  = MacPrmBr(B123U_Iv,IB1,R,dRdchi,dZdchi,T,IIT);
BIpF  = MacPrmBr(B123U_Ip,IB1,R,dRdchi,dZdchi,T,IIT);

BIpnF = MPRM * BIvF;

BIvR  = Fexp*BIvF;
BIpR  = Fexp*BIpF;
BIpnR = Fexp*BIpnF;

hf=figure(40);
CIv = V' * BIvF;  
plot(abs(CIv),'r-o','LineWidth',2,'MarkerSize',7), hold on,
xlabel('mode number','FontSize',16,'FontWeight','Bold')
ylabel('mode amplitude','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf=figure(41);
plot(tt*180/pi,abs(BIvR),'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(BIpnR),'r-','LineWidth',2), hold on,
%plot(tt,abs(BIpR),'r--','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('|b_n(n=1)| [Gauss/kA-turn]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('vacuum field','total response field')

d = -CEv(1)/CIv(1);
SOL_I = log10(sqrt(sum(abs(d*BIpnR+BEpnR).^2)*ht));
d = [abs(d) angle(d)*180/pi+120 d/abs(d) SOL_I];
SOLUTION_I = d 
hf=figure(42);
plot(tt*180/pi,abs(d(1)*BIvR*d(3)+BEvR),'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(d(1)*BIpnR*d(3)+BEpnR),'r-','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('|b_n(n=1)| [Gauss]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('vacuum field','total response field')

A       = linspace(0,1,101);
P       = linspace(0,360,37);
[AA,PP] = meshgrid(A,P);
PPD     = exp(i*(PP-120)*pi/180);
SOL     = AA;
for k=1:size(AA,1)
    for l=1:size(AA,2)
        SOL(k,l) = sum(abs(AA(k,l)*BIpnR*PPD(k,l)+BEpnR).^2);
    end
end
SOL = log10(sqrt(SOL*ht));

[Y,I] = min(SOL,[],2); [Y,K] = min(Y);
d = [AA(K,I(K)) PP(K,I(K)) PPD(K,I(K)) SOL(K,I(K))]; 
OPT_SOLUTION_I = d

hf=figure(43);
pcolor(AA,PP,SOL), hold on, shading interp
contour(AA,PP,SOL,20,'b-'), hold on,
xlabel('I_I[kA-turn]','FontSize',16,'FontWeight','Bold')
ylabel('I-coil Phase [degree]','FontSize',16,'FontWeight','Bold')
colorbar,  colormap(hot)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf=figure(44);
plot(tt*180/pi,abs(d(1)*BIvR*d(3)+BEvR),'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(d(1)*BIpnR*d(3)+BEpnR),'r-','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('|b_n(n=1)| [Gauss]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('vacuum field','total response field')

% combine both I-coil and C-coil
% an outer loop to do phase scan
PC = linspace(0,360,37);
PI = linspace(0,360,37);
PC = 60;
PI = 180;
[PPC,PPI] = meshgrid(PC,PI);
SOLP = PPC;

for kp=1:size(PPC,1)
for lp=1:size(PPC,2)
iteration = [kp lp]

% scan both I- and C-coils current (in kA-turn) 
% to try to cancel plasma response due to error field
% first define a criterion for the total response
% as int|BpR|dtheta
AC = linspace(-2,2,101);
AI = linspace(0,1,101);

[AAC,AAI] = meshgrid(AC,AI);
SOLA = AAC;

PC = exp(i*PPC(kp,lp)/180*pi);
PI = exp(i*(PPI(kp,lp)-120)/180*pi);

for k=1:size(AAC,1)
    for l=1:size(AAC,2)
        SOLA(k,l) = sum(abs(AAC(k,l)*BCpnR*PC + AAI(k,l)*BIpnR*PI + BEpnR).^2);
    end
end
SOLA = log10(sqrt(SOLA*ht));

[Y,I] = min(SOLA,[],2); [Y,K] = min(Y);
d = [AAC(K,I(K)) AAI(K,I(K)) SOLA(K,I(K))]; 
if size(PPC,1)==1, OPT_SOLUTION_A = d, end
SOLP(kp,lp) = d(3);
end
end

[Y,I] = min(SOLP,[],2); [Y,K] = min(Y);
d = [PPC(K,I(K)) PPI(K,I(K)) SOLP(K,I(K))]; 
OPT_SOLUTION_P = d

if size(PPC,1) == 1
hf=figure(50);
pcolor(AAC,AAI,SOLA), hold on, shading interp
contour(AAC,AAI,SOLA,20,'b-'), hold on,
xlabel('I_C[kA-turn]','FontSize',16,'FontWeight','Bold')
ylabel('I_I[kA-turn]','FontSize',16,'FontWeight','Bold')
colorbar,  colormap(hot)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf=figure(51);
d = OPT_SOLUTION_A;
plot(tt*180/pi,abs(d(1)*BCvR*PC+d(2)*BIvR*PI+BEvR),'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(d(1)*BCpnR*PC+d(2)*BIpnR*PI+BEpnR),'r-','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('|b_n(n=1)| [Gauss]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('vacuum field','total response field')

else 
hf=figure(52);
pcolor(PPC,PPI,SOLP), hold on, shading interp
contour(PPC,PPI,SOLP,20,'b-'), hold on,
xlabel('C-coil phase [degree]','FontSize',16,'FontWeight','Bold')
ylabel('Upper I-coil phase [degree]','FontSize',16,'FontWeight','Bold')
colorbar,  colormap(hot)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
end

% compare various solutions in complex plane
hf=figure(11);
plot(real(CEv(1)),imag(CEv(1)),'ro','LineWidth',3,'MarkerSize',9), hold on,
plot(real(CCv(1)),imag(CCv(1)),'r^','LineWidth',3,'MarkerSize',9), hold on,
plot(real(CIv(1)),imag(CIv(1)),'rv','LineWidth',3,'MarkerSize',9), hold on,
d = OPT_SOLUTION_C; c = -CCv(1)*d(1)*d(3);
plot(real(c),imag(c),'b+','LineWidth',3,'MarkerSize',9), hold on,
d = OPT_SOLUTION_I; c = -CIv(1)*d(1)*d(3);
plot(real(c),imag(c),'bx','LineWidth',3,'MarkerSize',9), hold on,
if size(PPC,1)==1
   d = OPT_SOLUTION_A; c = -(CCv(1)*d(1)*PC+CIv(1)*d(2)*PI);
   plot(real(c),imag(c),'bs','LineWidth',3,'MarkerSize',9), hold on,
end
xlabel('real axis','FontSize',16,'FontWeight','Bold')
ylabel('imaginary axis','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
if size(PPC,1)==1
   legend('EF','C-coil','I-coil','C-opt','I-opt','CI-opt')
else
   legend('EF','C-coil','I-coil','C-opt','I-opt')
end
a=axis;
plot([0 0],[a(3) a(4)],'k--'), hold on,
plot([a(1) a(2)],[0 0],'k--'), hold on,


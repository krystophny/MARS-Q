% generate and test plasma response matrices

global Mac

B0EXP= 0.4971;
mm   = [-29:29];
tt   = linspace(-pi,pi,201)';
MGG  = length(mm);

IB1 = 300;   %always associated with data B123U_*
IB2 = 127;   %always associated with data B123U2_*

Fexp = exp(i*tt*mm);
ht   = tt(2)-tt(1);

MacMainMAST_ALL_G
Mac.Mm = mm';
[T,IIT]=MacPrmGetT(R,Z,IB1);

SDIR = '/.automount/funsrv1/root/home/yliu/MAST_ADD/MAST26128/';
B123U_ALLv = load([SDIR 'B123U_ALLv']); 
B123U_ALLp = load([SDIR 'B123U_ALLp']); 
B123U_Cv   = load([SDIR 'B123U_Cv']); 
B123U_Cp   = load([SDIR 'B123U_Cp']); 
B123U_Ev   = load([SDIR 'B123U_Ev']); 
B123U_Ep   = load([SDIR 'B123U_Ep']); 

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
plot(diag(S),'ro','LineWidth',2,'MarkerSize',7), hold on,
xlabel('mode number','FontSize',16,'FontWeight','Bold')
ylabel('SVD singular value','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

% compute the radial component of the vacuum n=1 error field
% and test the true error field response
BEvF  = MacPrmBr(B123U_Ev,IB1,R,dRdchi,dZdchi,T,IIT);
BEpF  = MacPrmBr(B123U_Ep,IB1,R,dRdchi,dZdchi,T,IIT);

BEpnF = MPRM * BEvF;

BEvR  = Fexp*BEvF;
BEpR  = Fexp*BEpF;
BEpnR = Fexp*BEpnF;

hf=figure(20);
CEv = V' * BEvF;  
plot(abs(CEv),'r-o','LineWidth',2,'MarkerSize',7), hold on,
xlabel('mode number','FontSize',16,'FontWeight','Bold')
ylabel('mode amplitude (EF)','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf=figure(21);
plot(tt*180/pi,abs(BEvR)*B0EXP,'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(BEpR)*B0EXP,'r-','LineWidth',2), hold on,
%plot(tt*180/pi,abs(BEpnR),'r-','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('EF |b_n(n=1)| (Gauss)','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
%legend('b_n^v','b_n^p','M*b_n^v')
legend('b_n^v','b_n^p')

% test the EFCC response
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
ylabel('mode amplitude (EFCC)','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf=figure(31);
plot(tt*180/pi,abs(BCvR),'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(BCpR),'r--','LineWidth',2), hold on,
plot(tt*180/pi,abs(BCpnR),'r-','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('|b_n(n=1)| (EFCC)','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('b_n^v','b_n^p','M*b_n^v')

%optimal correction by correcting the most sensitive SVD component of vacuum field
d = -CEv(1)/CCv(1);
SOL_C = log10(sqrt(sum(abs(d*BCpR+BEpR).^2)*ht));
dsol = [real(d) imag(d) SOL_C];
OPT_SOLUTION_D = dsol     %optimal solution based on SVD criterion 
d_SVD = d;
hf=figure(32);
plot(tt*180/pi,abs(d*BCvR+BEvR)*B0EXP,'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(d*BCpR+BEpR)*B0EXP,'r-','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('EF+EFCC |b_n(n=1)| (Gauss)','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('vacuum','response')

%find optimal solution by minimizing ||d*BCpR+BEpR|| w.r.t. d ==> criterion E
A       = linspace(-2e-3,2e-3,51);
P       = A;
[AA,PP] = meshgrid(A,P);
dd      = AA + PP*i;
SOL     = AA;
for k=1:size(AA,1)
    for l=1:size(AA,2)
        SOL(k,l) = sum(abs(dd(k,l)*BCpR+BEpR).^2);
    end
end
SOL = log10(sqrt(SOL*ht));

[Y,I] = min(SOL,[],2); [Y,K] = min(Y);  d = dd(K,I(K));
dsol = [d SOL(K,I(K))]; 
OPT_SOLUTION_E = dsol     %optimal solution based on full field correction
d_FULL = d;

hf=figure(33);
pcolor(AA,PP,SOL), hold on, shading interp
contour(AA,PP,SOL,20,'b-'), hold on,
xlabel('Re(I_{EFCC})','FontSize',16,'FontWeight','Bold')
ylabel('Im(I_{EFCC})','FontSize',16,'FontWeight','Bold')
colorbar,  colormap(hot)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf=figure(34);
plot(tt*180/pi,abs(d*BCvR+BEvR),'b--','LineWidth',2), hold on,
plot(tt*180/pi,abs(d*BCpR+BEpR),'r-','LineWidth',2), hold on,
xlabel('geometric poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
ylabel('|b_n(n=1)| (FULL:EF+EFCC)','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('vacuum','response')

%convert optimal EFCC current to MAST coil current
  deltaP = 83/180*pi;
  R0EXP  = 0.7712;
  B0EXP  = 0.4971;
  mu0    = 4e-7*pi;
  n      = -1;
  facI0  = R0EXP*B0EXP/mu0*1e-3;
  facI1  = sqrt(2)*sin(n*deltaP/2)/n/pi;
  facI2  =-1.0;

  In     = d_SVD;
  Ic1SVD = (imag(In/facI2)+real(In/facI2))/2*facI0/facI1; %[kAt]  
  Ic2SVD = (imag(In/facI2)-real(In/facI2))/2*facI0/facI1; %[kAt]

  In     = d_FULL;
  Ic1FULL= (imag(In/facI2)+real(In/facI2))/2*facI0/facI1; %[kAt]  
  Ic2FULL= (imag(In/facI2)-real(In/facI2))/2*facI0/facI1; %[kAt]

  IEFCC = [Ic1SVD Ic2SVD Ic1FULL Ic2FULL]




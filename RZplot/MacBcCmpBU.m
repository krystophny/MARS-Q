%%Compare B2U/B1U as a function of gamma

format short e
global Fc A1 A2a A2b A3

kcase = 4;

% CARIDDIcirc from CARIDDI
if kcase == 2
  Fc.DIR0    = '/.automount/funsrv1/root/home/yliu/DataMarsf/Backward/2dWall/';
  Fc.cariddi = 1;   
  Fc.kpsi    =-1;
  Fc.tauA = 1.0587e-6;
  Fc.mm  = -5:9;
  Fc.rs  = 1.09665;
  Fc.CC  = 0;  %include control coils ?
  Fc.fac   = Mac.R0EXP/(4e-7*pi);
  Fc.ks  = [1 5 8 4 1 8 4];
  kk  = -5:9;
  Fc.SS    = 'b-s';
  Fc.DIR   = [Fc.DIR0 'Passive/'];
  MacBcCar2Mar([Fc.DIR0 'quintuples_eqv_moreholes.mat'])
  gamma2p  = Fc.gamma2*Fc.tauA;
  if Fc.CC > 0
     psisp     = Fc.psis;
     Fc.DIR   = [Fc.DIR0 'Active/'];
     MacBcCar2Mar([Fc.DIR0 'quintuples_f_moreholes.mat'])
     gamma2a = Fc.gamma2*Fc.tauA;
     psisa   = Fc.psis*Fc.fac;
  end
end

% ITER from CARIDDI
if kcase == 4
  Fc.DIR0    = '/.automount/funsrv1/root/home/yliu/Scen4_V15/Backward/';
  Fc.cariddi = 1;   
  Fc.kpsi    =-1;
  Fc.tauA = 0.6454E-06;
  Fc.mm  = -9:39;
  Fc.rs  = 1.01;
  Fc.CC  = 0;
  Fc.fac   = Mac.R0EXP/(4e-7*pi);
  kk  = -9:39;
  if Fc.cariddi==1
     Fc.SS    = 'b-s';
     Fc.ks  = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
     Fc.DIR   = [Fc.DIR0 'Passive/'];
     MacBcCar2Mar([Fc.DIR0 'quadruples_ITER_new_3.mat'])
  else
     Fc.SS    = 'r-o';
     Fc.ks  = [1 2 3 4 5 6];
     Fc.DIR   = [Fc.DIR0 'Marsf/'];
     Fc.gamma2 = [6.454e-8 6.454e-7 6.454e-6 6.454e-5 6.454e-4 6.454e-3 6.454e-1]/Fc.tauA;
  end
  gamma2p  = Fc.gamma2*Fc.tauA;
  if Fc.CC > 0
     psisp     = Fc.psis;
     Fc.DIR   = [Fc.DIR0 'Active/'];
     MacBcCar2Mar([Fc.DIR0 'quintuples_f_moreholes.mat'])
     gamma2a = Fc.gamma2*Fc.tauA;
     psisa   = Fc.psis*Fc.fac;
  end
end

Fc.ii = kk - Fc.mm(1) + 1;
Fc.plotA = 3;  %1,2,3
Fc.plotB = 5;  %4,5
Fc.plotD = 0;
Fc.plotC = 0;
Fc.plotS = 6;

res_gammap = [[1:length(gamma2p)]' gamma2p(:)] 

%transform CARIDDI output format to our pre-defined data format 
%coordinate transform matrices from CARIDDI geometrical theta to MARS-F chi
if Fc.cariddi > 0
  [A1]=MacBcGetCtMatN3(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
  [A2a,A2b,A3]=MacBcGetCtMatT3(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
  %N=size(A1,1); A1=eye(N); A2a=A1; A2b=zeros(N,N); A3=A1;
  %A1 = A1*1.25;
end
  
%read input matrices
DIR = Fc.DIR;

BU21 = zeros(length(kk),length(gamma2p));
for jj=1:length(gamma2p)
    s0 = gamma2p(jj);
    g0 = num2str(s0);
    [B1U0,B2U0,B3U0]=MacBcReadMat([DIR 'BU' g0 '_ALL'],1);
    V = diag(B2U0);
    BU21(:,jj) = V(:);
end


if Fc.plotB == 5
m0=1;
figure(m0+10*Fc.plotB)
k0=m0 - Fc.mm(1) + 1;
%semilogx(gamma2p,real(BU21(k0,:)),Fc.SS,'LineWidth',2), hold on,
semilogx(gamma2p,imag(BU21(k0,:)),Fc.SS,'LineWidth',0.5), hold on,
xlabel('\gamma_0'), ylabel('B2U/B1U')

m0=2;
figure(m0+10*Fc.plotB)
k0=m0 - Fc.mm(1) + 1;
%semilogx(gamma2p,real(BU21(k0,:)),Fc.SS,'LineWidth',2), hold on,
semilogx(gamma2p,imag(BU21(k0,:)),Fc.SS,'LineWidth',0.5), hold on,
xlabel('\gamma_0'), ylabel('B2U/B1U')

m0=3;
figure(m0+10*Fc.plotB)
k0=m0 - Fc.mm(1) + 1;
%semilogx(gamma2p,real(BU21(k0,:)),Fc.SS,'LineWidth',2), hold on,
semilogx(gamma2p,imag(BU21(k0,:)),Fc.SS,'LineWidth',0.5), hold on,
xlabel('\gamma_0'), ylabel('B2U/B1U')

m0=4;
figure(m0+10*Fc.plotB)
k0=m0 - Fc.mm(1) + 1;
%semilogx(gamma2p,real(BU21(k0,:)),Fc.SS,'LineWidth',2), hold on,
semilogx(gamma2p,imag(BU21(k0,:)),Fc.SS,'LineWidth',0.5), hold on,
xlabel('\gamma_0'), ylabel('B2U/B1U')

m0=5;
figure(m0+10*Fc.plotB)
k0=m0 - Fc.mm(1) + 1;
%semilogx(gamma2p,real(BU21(k0,:)),Fc.SS,'LineWidth',2), hold on,
semilogx(gamma2p,imag(BU21(k0,:)),Fc.SS,'LineWidth',0.5), hold on,
xlabel('\gamma_0'), ylabel('B2U/B1U')
end


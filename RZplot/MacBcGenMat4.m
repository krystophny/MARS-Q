%function MacBcGenMat(R,Z,jacobian)
%%Generate input coupling matrices for backward coupling scheme
%%Input should be created by CARIDDI
%%with active coils (Fc.CC>0)
%%need to run MacMain*.m first to generate (R,Z,jacobian) 
%%using triple-data with gamma0=0 and arbitrary (gamma1,gamma2)

format short e
global Fc A1 A2a A2b A3

kcase = 8;

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
  kk  = -3:15;
  if Fc.cariddi==1
     Fc.SS    = 'b-+';
     Fc.ks  = [1 9 11 14];
     Fc.DIR   = [Fc.DIR0 'Passive/'];
     MacBcCar2Mar([Fc.DIR0 'quadruples_ITER_new_3.mat'])
  else
     Fc.SS    = 'r-o';
     Fc.ks  = [3 4 5 7];
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

% JET68875 data from CARIDDI
if kcase == 8
  Fc.DIR0    = '/.automount/funsrv1/root/home/yliu/JET68875_CARMA/';
  Fc.cariddi = 1;   
  Fc.kpsi    =-1;
  Fc.tauA = 5.0390e-7;
  Fc.mm  = -9:25;
  Fc.rs  = 1.05;
  Fc.CC  = 0;  %include control coils ?
  Fc.fac   = Mac.R0EXP/(4e-7*pi);
  Fc.ks  = [2 8 13 8 1 8 4];
  kk  = -9:25;
  Fc.SS    = 'b-s';
  Fc.DIR   = [Fc.DIR0 'Passive/'];
  MacBcCar2Mar([Fc.DIR0 'quadruples01.mat'])
  gamma2p  = Fc.gamma2*Fc.tauA;
  if Fc.CC > 0
     psisp     = Fc.psis;
     Fc.DIR   = [Fc.DIR0 'Active/'];
     MacBcCar2Mar([Fc.DIR0 'quadruples05.mat'])
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
s0 = gamma2p(Fc.ks(1));
s1 = gamma2p(Fc.ks(2));
s2 = gamma2p(Fc.ks(3));
s3 = gamma2p(Fc.ks(4));  %for test
g0 = num2str(s0);
g1 = num2str(s1);
g2 = num2str(s2);
g3 = num2str(s3);

%transform CARIDDI output format to our pre-defined data format 
%coordinate transform matrices from CARIDDI geometrical theta to MARS-F chi
if Fc.cariddi > 0
  [A1]=MacBcGetCtMatN3(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
  [A2a,A2b,A3]=MacBcGetCtMatT3(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
end
  
%read input matrices
DIR = Fc.DIR;
[B1U0,B2U0,B3U0,B2U0A,B3U0A]=MacBcReadMat([DIR 'BU' g0 '_ALL'],1);
[B1U1,B2U1,B3U1,B2U1A,B3U1A]=MacBcReadMat([DIR 'BU' g1 '_ALL'],1);
[B1U2,B2U2,B3U2,B2U2A,B3U2A]=MacBcReadMat([DIR 'BU' g2 '_ALL'],1);
[B1U3,B2U3,B3U3]=MacBcReadMat([DIR 'BU' g3 '_ALL'],1);

N = length(Fc.ii);
A20 = eye(N);
B20 = B2U0;
A21 = ((B2U1A-B2U0A)/s1-(B2U2A-B2U0A)/s2)*inv(B2U2A-B2U1A);
B21 = A21*B2U1 + (B2U1-B2U0)/s1;

%compose a big matrix system for toroidal component
A = [ s0*B3U0  s1*B3U1  s2*B3U2
        -B1U0    -B1U1    -B1U2
     -s0*B1U0 -s1*B1U1 -s2*B1U2];
b = -[B3U0 B3U1 B3U2];

%add active coil part
if Fc.CC > 0
   A = [A [s4*B3U4 s5*B3U5; -B1U4 -B1U5; -s4*B1U4 -s5*B1U5]
        MA0 MA0 MA0 MA1 MA1
        MA0 MA0 MA0 s4*MA1 s5*MA1];
   b = [b -B3U4 -B3U5];
end

Ainv = pinv(A);
u = b*Ainv;
N = length(Fc.ii);
A30 = eye(N);
A31 = u(:,1:N);
B30 = u(:,N+1:2*N);
B31 = u(:,2*N+1:3*N);
if Fc.CC > 0
   C30 = u(:,3*N+1:3*N+MK);
   C31 = u(:,3*N+MK+1:end);
end

%sensor matrices
if Fc.CC > 0
   u    = bs*Ainv;
   ML   = size(u,1);
   A30S = zeros(ML,N);
   A31S = u(:,1:N);
   B30S = u(:,N+1:2*N);
   B31S = u(:,2*N+1:3*N);
   C30S = u(:,3*N+1:3*N+MK);
   C31S = u(:,3*N+MK+1:end);
end

% test coupling matrices
condnum_tor = cond(A);

%response matrices for sensors

res_condnum = [condnum_pol condnum_tor]


k = 1 - Fc.mm(1) + 1; %m=1 diagonal element
res_field = [s0 B2U0(k,k)/B1U0(k,k) B3U0(k,k)/B1U0(k,k) 
             s1 B2U1(k,k)/B1U1(k,k) B3U1(k,k)/B1U1(k,k) 
             s2 B2U2(k,k)/B1U2(k,k) B3U2(k,k)/B1U2(k,k) 
             s3 B2U3(k,k)/B1U3(k,k) B3U3(k,k)/B1U3(k,k)]

D2e = ( (A20+s3*A21)*B2U3 - (B20+s3*B21)*B1U3 )/max(max(abs(B1U3)));
D3e = ( (A30+s3*A31)*B3U3 - (B30+s3*B31)*B1U3 )/max(max(abs(B1U3)));
res = [max(max(abs(D2e))) max(max(abs(D3e)))];

if Fc.CC > 0
   D2c = ( (A20+s6*A21)*B2U6 - (B20+s6*B21)*B1U6 + (C20+s6*C21) )./abs(B1U6);
   D3c = ( (A30+s6*A31)*B3U6 - (B30+s6*B31)*B1U6 + (C30+s6*C31) )./abs(B1U6);

   S2e = ( (A20S+s3*A21S)*B2U3 - (B20S+s3*B21S)*B1U3 - psisp(:,:,Fc.ks(4)) )./abs(psisp(:,:,Fc.ks(4)));
   S3e = ( (A30S+s3*A31S)*B3U3 - (B30S+s3*B31S)*B1U3 - psisp(:,:,Fc.ks(4)) )./abs(psisp(:,:,Fc.ks(4)));

   S2c = ( (A20S+s6*A21S)*B2U6 - (B20S+s6*B21S)*B1U6 + (C20S+s6*C21S) - psisa(:,:,Fc.ks(7)) )./abs(psisa(:,:,Fc.ks(7)));
   S3c = ( (A30S+s6*A31S)*B3U6 - (B30S+s6*B31S)*B1U6 + (C30S+s6*C31S) - psisa(:,:,Fc.ks(7)) )./abs(psisa(:,:,Fc.ks(7)));

   res = [res
          max(max(abs(D2c))) max(max(abs(D3c)))
          max(max(abs(S2e))) max(max(abs(S3e))) 
          max(max(abs(S2c))) max(max(abs(S3c)))];
end

res_test = res

%check maximum values of coupling matrices
res = [max(max(abs(A20))) max(max(abs(A21)))
       max(max(abs(A30))) max(max(abs(A31)))
       max(max(abs(B20))) max(max(abs(B21)))
       max(max(abs(B30))) max(max(abs(B31)))];
if Fc.CC > 0
res = [res [max(max(abs(A20S))) max(max(abs(A21S)));...
            max(max(abs(A30S))) max(max(abs(A31S)));...
            max(max(abs(B20S))) max(max(abs(B21S)));...
            max(max(abs(B30S))) max(max(abs(B31S)));]
       max(max(abs(C20))) max(max(abs(C21))) max(max(abs(C20S))) max(max(abs(C21S))) 
       max(max(abs(C30))) max(max(abs(C31))) max(max(abs(C30S))) max(max(abs(C31S)))];
end
res_matmax = res 
   
%save coupling matrices
ress = [real(A20) imag(A20) real(A21) imag(A21)...
        real(B20) imag(B20) real(B21) imag(B21)...
        real(A30) imag(A30) real(A31) imag(A31)...
        real(B30) imag(B30) real(B31) imag(B31)];
if Fc.CC > 0
   ress = [ress...
           real(C20) imag(C20) real(C21) imag(C21)...
           real(C30) imag(C30) real(C31) imag(C31)]; 
end
     
save BWCMAT ress -ascii -double
eval(['!mv BWCMAT ' Fc.DIR0])

%save sensor matrices
if Fc.CC > 0
   ress = [real(A20S) imag(A20S) real(A21S) imag(A21S)...
           real(B20S) imag(B20S) real(B21S) imag(B21S)...
           real(A30S) imag(A30S) real(A31S) imag(A31S)...
           real(B30S) imag(B30S) real(B31S) imag(B31S)...
           real(C20S) imag(C20S) real(C21S) imag(C21S)...
           real(C30S) imag(C30S) real(C31S) imag(C31S)]; 

   save BWCMATS ress -ascii -double
   eval(['!mv BWCMATS ' Fc.DIR0])
end
     
% plotting...
if Fc.plotA == 1
figure(1+10*Fc.plotA)
plot(kk,real(B20),Fc.SS), hold on,
xlabel('m'), ylabel('Re(B20)')

figure(2+10*Fc.plotA)
plot(kk,real(A21),Fc.SS), hold on,
xlabel('m'), ylabel('Re(A21)')

figure(3+10*Fc.plotA)
plot(kk,real(B21),Fc.SS), hold on,
xlabel('m'), ylabel('Re(B21)')

figure(4+10*Fc.plotA)
plot(kk,real(B30),Fc.SS), hold on,
xlabel('m'), ylabel('Re(B30)')

figure(5+10*Fc.plotA)
plot(kk,real(A31),Fc.SS), hold on,
xlabel('m'), ylabel('Re(A31)')

figure(6+10*Fc.plotA)
plot(kk,real(B31),Fc.SS), hold on,
xlabel('m'), ylabel('Re(B31)')
end

if Fc.plotA == 2
figure(1+10*Fc.plotA)
plot(kk,imag(B20),Fc.SS), hold on,
xlabel('m'), ylabel('Im(B20)')

figure(2+10*Fc.plotA)
plot(kk,imag(A21),Fc.SS), hold on,
xlabel('m'), ylabel('Im(A21)')

figure(3+10*Fc.plotA)
plot(kk,imag(B21),Fc.SS), hold on,
xlabel('m'), ylabel('Im(B21)')

figure(4+10*Fc.plotA)
plot(kk,imag(B30),Fc.SS), hold on,
xlabel('m'), ylabel('Im(B30)')

figure(5+10*Fc.plotA)
plot(kk,imag(A31),Fc.SS), hold on,
xlabel('m'), ylabel('Im(A31)')

figure(6+10*Fc.plotA)
plot(kk,imag(B31),Fc.SS), hold on,
xlabel('m'), ylabel('Im(B31)')
end

% plotting...
if Fc.plotA == 3
figure(1+10*Fc.plotA)
subplot(3,2,1)
plot(kk,abs(diag(B20)),Fc.SS,'MarkerSize',7), hold on,
xlabel('m'), ylabel('diag(|B20|)')

subplot(3,2,2)
plot(kk,abs(diag(B30)),Fc.SS,'MarkerSize',7), hold on,
xlabel('m'), ylabel('diag(|B30|)')

subplot(3,2,3)
plot(kk,abs(diag(A21)),Fc.SS,'MarkerSize',7), hold on,
xlabel('m'), ylabel('diag(|A21|)')

subplot(3,2,4)
plot(kk,abs(diag(A31)),Fc.SS,'MarkerSize',7), hold on,
xlabel('m'), ylabel('diag(|A31|)')

subplot(3,2,5)
plot(kk,abs(diag(B21)),Fc.SS,'MarkerSize',7), hold on,
xlabel('m'), ylabel('diag(|B21|)')

subplot(3,2,6)
plot(kk,abs(diag(B31)),Fc.SS,'MarkerSize',7), hold on,
xlabel('m'), ylabel('diag(|B31|)')
end

if Fc.plotB == 4
figure(1+10*Fc.plotB)
semilogy(kk,abs(diag(B1U1)),Fc.SS), hold on,
xlabel('m'), ylabel('|B1U1|')

figure(2+10*Fc.plotB)
semilogy(kk,abs(diag(B2U1)),Fc.SS), hold on,
xlabel('m'), ylabel('|B2U1|')

figure(3+10*Fc.plotB)
semilogy(kk,abs(diag(B3U1)),Fc.SS), hold on,
xlabel('m'), ylabel('|B3U1|')
end

if Fc.plotB == 5
figure(1+10*Fc.plotB)
plot(kk,real(diag(B2U1)),Fc.SS,'LineWidth',2), hold on,
plot(kk,imag(diag(B2U1)),Fc.SS), hold on,
xlabel('m'), ylabel('B2U1')

figure(2+10*Fc.plotB)
plot(kk,real(diag(B2U2)),Fc.SS,'LineWidth',2), hold on,
plot(kk,imag(diag(B2U2)),Fc.SS), hold on,
xlabel('m'), ylabel('B2U2')

figure(3+10*Fc.plotB)
plot(kk,real(diag(B2U3)),Fc.SS,'LineWidth',2), hold on,
plot(kk,imag(diag(B2U3)),Fc.SS), hold on,
xlabel('m'), ylabel('B2U3')
end

%check vacuum response of coils to sensors
if Fc.CC > 0
   kc0  = 1;  %# of active coil
   psis = zeros(size(psisa,1),size(psisa,3));
   for k=1:size(psisa,3)
       psis(:,k) = psisa(:,kc0,k)/Fc.fac;
   end
   if Fc.plotS > 0
      N = round((size(psisa,3)-1)/2);
      figure(0 + 10*Fc.plotS)
      semilogx(abs(gamma2a(2:N+1))/Fc.tauA,real(psis(1,2:N+1)),'r-o','LineWidth',1,'MarkerSize',8), hold on,
      semilogx(abs(gamma2a(N+2:end))/Fc.tauA,real(psis(1,N+2:end)),'r--s','LineWidth',1,'MarkerSize',8), hold on,

      figure(1 + 10*Fc.plotS)
      semilogx(abs(gamma2a(2:N+1))/Fc.tauA,imag(psis(1,2:N+1)),'r-o','LineWidth',1,'MarkerSize',8), hold on,
      semilogx(abs(gamma2a(N+2:end))/Fc.tauA,imag(psis(1,N+2:end)),'r--s','LineWidth',1,'MarkerSize',8), hold on,

		psis1 = [psis(1,1) psis(1,N+2:end)];
      psis1 = psis1/psis1(1);
      figure(1)
      plot(real(psis1),imag(psis1),'k-'), hold on,
   end
   psis0 = psis(1,1)*Fc.fac
end


%function MacBcGenMat(R,Z,jacobian)
%%Generate input coupling matrices for backward coupling scheme
%%Input should be created by CARIDDI
%%no active coils
%%need to run MacMain*.m first to generate (R,Z,jacobian) 
%%using triple-data with arbitrary (gamma1,gamma2,gamma3)

format short e
global Fc A1 A2a A2b A3

kcase = 8;

% CARIDDIcirc from CARIDDI
if kcase == 2
  Fc.DIR = '/.automount/funsrv1/root/home/yliu/CARIDDIcirc/DataCariddi/';
  Fc.cariddi = 1;   
  Fc.tauA = 1.0587e-6;
  Fc.mm  = -5:9;
  Fc.rs  = 1.09665;
  kk  = -5:9;
  Fc.SS    = 'b-s';
  MacBcCar2Mar([Fc.DIR 'quadruples3_moreholes.mat'])
end

% ITER from CARIDDI
if kcase == 4
  Fc.DIR = '/.automount/funsrv1/root/home/yliu/Scen4_V15/Backward/';
  Fc.cariddi = 1;
  Fc.kpsi    = -1; %include sensor signals for feedback case   
  Fc.tauA = 0.6454E+2;
  Fc.mm  = -9:39;
  Fc.rs  = 1.01;
  kk  = -9:20;
  Fc.SS    = 'b-o';
  MacBcCar2Mar([Fc.DIR 'quadruples_ITER_new_1.mat'])
end

% JET from CARIDDI
if kcase == 8
  Fc.DIR = '/.automount/funsrv1/root/home/yliu/JET68875_CARMA/';
  Fc.cariddi = 1;
  Fc.kpsi    = -1; %include sensor signals for feedback case   
  Fc.tauA = 5.0390e-07;
  Fc.mm  = -9:25;
  Fc.rs  = 1.05;
  kk  = -9:25;
  Fc.SS    = 'k-s';
  if Fc.cariddi>0, MacBcCar2Mar([Fc.DIR 'quadruples01.mat']), end
end

Fc.ii = kk - Fc.mm(1) + 1;
Fc.plotA = 3;  %1,2,3
Fc.plotB = 5;  %4,5
Fc.plotD = 0;
Fc.plotC = 0;

DIR = Fc.DIR;

if Fc.cariddi >0
gamma2 = Fc.gamma2*Fc.tauA; %normalised by tauA
res_gamma = [[1:length(gamma2)]' gamma2(:)] 
s0 = gamma2(2);
s1 = gamma2(8 );
s2 = gamma2(13);
s3 = gamma2(8);
g0 = num2str(s0);
g1 = num2str(s1);
g2 = num2str(s2);
g3 = num2str(s3);
else
gamma2 = [5.0390e-7 1.5117e-4 5.0390e-1];
s0 = gamma2(1);
s1 = gamma2(2 );
s2 = gamma2(3);
s3 = gamma2(2);
g0 = num2str(s0);
g1 = num2str(s1);
g2 = num2str(s2);
g3 = num2str(s3);
end

%transform CARIDDI output format to our pre-defined data format 
%coordinate transform matrices from CARIDDI geometrical theta to MARS-F chi
if Fc.cariddi > 0
  [A1]=MacBcGetCtMatN3(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
  [A2a,A2b,A3]=MacBcGetCtMatT3(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
  %N=size(A1,1); A1=eye(N); A2a=A1; A2b=zeros(N,N); A3=A1;
  %A1 = A1*1.25;
end
  
%read input matrices
if Fc.cariddi>0
[B1U0,B2U0,B3U0]=MacBcReadMat([DIR 'BU' g0 '_ALL'],1);
[B1U1,B2U1,B3U1]=MacBcReadMat([DIR 'BU' g1 '_ALL'],1);
[B1U2,B2U2,B3U2]=MacBcReadMat([DIR 'BU' g2 '_ALL'],1);
else
[B1U0,B2U0,B3U0]=MacBcReadMat([DIR 'MarsF/BU' g0 '_ALL'],1);
[B1U1,B2U1,B3U1]=MacBcReadMat([DIR 'MarsF/BU' g1 '_ALL'],1);
[B1U2,B2U2,B3U2]=MacBcReadMat([DIR 'MarsF/BU' g2 '_ALL'],1);
end

%normalise field harmonics
D0 = kk; I = find(D0==0); if length(I)==1, D0(I)=1; end
D0 = diag(D0);   sn=1/s1;
D1 = D0/6;       
D2 = abs(D0)/8;   
D3 = D0/0.1; 

sn = 1.0; N = length(kk); D1=eye(N); D2=D1; D3=D1;    

%compose a big matrix system for poloidal component
A = [ sn*s0*D2*B2U0  sn*s1*D2*B2U1  sn*s2*D2*B2U2
           -D1*B1U0       -D1*B1U1       -D1*B1U2
     -sn*s0*D1*B1U0 -sn*s1*D1*B1U1 -sn*s2*D1*B1U2];
b = -[B2U0 B2U1 B2U2];
u = b*inv(A);
N = round(size(u,2)/3);
A20 = eye(N);
A21 = u(:,1:N)*D2*sn;
B20 = u(:,N+1:2*N)*D1;
B21 = u(:,2*N+1:end)*D1*sn;
% test coupling matrices
condnum_pol = cond(A);

%compose a big matrix system for toroidal component
A = [ sn*s0*D3*B3U0  sn*s1*D3*B3U1  sn*s2*D3*B3U2
           -D1*B1U0       -D1*B1U1       -D1*B1U2
     -sn*s0*D1*B1U0 -sn*s1*D1*B1U1 -sn*s2*D1*B1U2];
b = -[B3U0 B3U1 B3U2];
u = b*inv(A);
N = round(size(u,2)/3);
A30 = eye(N);
A31 = u(:,1:N)*D3*sn;
B30 = u(:,N+1:2*N)*D1;
B31 = u(:,2*N+1:end)*D1*sn;
% test coupling matrices
condnum_tor = cond(A);

res_condnum = [condnum_pol condnum_tor]

if Fc.cariddi>0
[B1U3,B2U3,B3U3]=MacBcReadMat([DIR 'BU' g3 '_ALL'],1);
else
[B1U3,B2U3,B3U3]=MacBcReadMat([DIR 'MarsF/BU' g3 '_ALL'],1);
end

k = 1 - Fc.mm(1) + 1; %m=1 diagonal element
res_field = [s0 B2U0(k,k)/B1U0(k,k) B3U0(k,k)/B1U0(k,k) 
             s1 B2U1(k,k)/B1U1(k,k) B3U1(k,k)/B1U1(k,k) 
             s2 B2U2(k,k)/B1U2(k,k) B3U2(k,k)/B1U2(k,k) 
             s3 B2U3(k,k)/B1U3(k,k) B3U3(k,k)/B1U3(k,k)]

D2e = ( (A20+s3*A21)*B2U3 - (B20+s3*B21)*B1U3 )./max(max(abs(B1U3)));
D3e = ( (A30+s3*A31)*B3U3 - (B30+s3*B31)*B1U3 )./max(max(abs(B1U3)));
res_test = [max(max(abs(D2e))) max(max(abs(D3e)))]

%save coupling matrices
ress = [real(A20) imag(A20) real(A21) imag(A21)...
        real(B20) imag(B20) real(B21) imag(B21)...
        real(A30) imag(A30) real(A31) imag(A31)...
        real(B30) imag(B30) real(B31) imag(B31)];
save BWCMAT ress -ascii -double
eval(['!mv BWCMAT ' DIR])

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
plot(kk,real(diag(B1U1*D1)),Fc.SS,'LineWidth',2), hold on,
plot(kk,imag(diag(B1U1*D1)),Fc.SS), hold on,
xlabel('m'), ylabel('B1U1')

figure(2+10*Fc.plotB)
plot(kk,real(diag(B2U1*D2)),Fc.SS,'LineWidth',2), hold on,
plot(kk,imag(diag(B2U1*D2)),Fc.SS), hold on,
xlabel('m'), ylabel('B2U1')

figure(3+10*Fc.plotB)
plot(kk,real(diag(B3U1*D3)),Fc.SS,'LineWidth',2), hold on,
plot(kk,imag(diag(B3U1*D3)),Fc.SS), hold on,
xlabel('m'), ylabel('B3U1')
end

if Fc.plotD>0
figure(1+10*Fc.plotD)
semilogy(kk,abs(D2),Fc.SS), hold on,
xlabel('m'), ylabel('|D2|')

figure(2+10*Fc.plotD)
semilogy(kk,abs(D3),Fc.SS), hold on,
xlabel('m'), ylabel('|D3|')
end


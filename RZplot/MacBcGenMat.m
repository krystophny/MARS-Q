%function MacBcGenMat(R,Z,jacobian)
%%Generate input coupling matrices for backward coupling scheme
%%need to run MacMain*.m first to generate (R,Z,jacobian) 

format short e
global Fc A1 A2a A2b A3 

kcase = 3;

% CARIDDIcirc from MARS-F
if kcase == 1
  Fc.DIR = '/.automount/funsrv1/root/home/yliu/CARIDDIcirc/DataMarsf/Backward/';
  Fc.cariddi = 0;   
  Fc.g2  = '2e-3i';
  Fc.g3  = '2e-4i';
  Fc.mm  = -5:9;
  Fc.rs  = 1.09665;
  kk  = Fc.mm;
  Fc.SS    = 'r-o';
end

% ITER from MARS-F
if kcase == 3
% Fc.DIR = '/.automount/funsrv1/root/home/yliu/ITER/Scen4_V14/DataMar/M2_20/';
  Fc.DIR = '/.automount/funsrv1/root/home/yliu/ITER/Scen4_V14/DataCar/';
  Fc.cariddi = 0;   
  Fc.g2  = '6.4540e-6';
  Fc.g3  = '6.4540e-6';
  Fc.mm  = -9:39;
  Fc.rs  = 1.01;
  kk  = Fc.mm;
  Fc.SS    = 'r-o';
end

Fc.ii = kk - Fc.mm(1) + 1;
Fc.plotA = 3; %1,2,3
Fc.plotB = 0; %4,5
Fc.plotD = 0;
Fc.plotC = 0;

DIR = Fc.DIR;
g2 = Fc.g2
g3 = Fc.g3
eval(['s2=' g2 ';']);
eval(['s3=' g3 ';']);

%transform CARIDDI output format to our pre-defined data format 
%coordinate transform matrices from CARIDDI geometrical theta to MARS-F chi
if Fc.cariddi > 0
  [A1,A2a,A2b,A3]=MacBcGetCtMat(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
end
  
sn = 1.0; N = length(kk); D1=eye(N); D2=D1; D3=D1;    

%read input matrices
[B1U0,B2U0,B3U0]=MacBcReadMat([DIR 'BU0_ALL']);
[B1U1,B2U1,B3U1]=MacBcReadMat([DIR 'BUinf_ALL']);
[B1U2,B2U2,B3U2]=MacBcReadMat([DIR 'BU' g2 '_ALL']);
 
%invert some matrices
B1U0v = inv(B1U0);
B1U1v = inv(B1U1);

%1st order Pade coupling matrices 
A20 = eye(size(B1U0));
B20 = A20*B2U0*B1U0v;
A21 = (B20*B1U2-B2U2)*inv(B2U2-B2U1*B1U1v*B1U2)/s2;
B21 = A21*B2U1*B1U1v;

A30 = eye(size(B1U0));
B30 = A30*B3U0*B1U0v;
A31 = (B30*B1U2-B3U2)*inv(B3U2-B3U1*B1U1v*B1U2)/s2;
B31 = A31*B3U1*B1U1v;

% test coupling matrices
condnum = [cond(B1U0) cond(B1U1) cond(B1U2)]

[B1U3,B2U3,B3U3]=MacBcReadMat([DIR 'BU' g3 '_ALL']);

D2e = ( (A20+s3*A21)*B2U3 - (B20+s3*B21)*B1U3 )./abs(B1U3);
D3e = ( (A30+s3*A31)*B3U3 - (B30+s3*B31)*B1U3 )./abs(B1U3);
res = [max(max(abs(D2e))) max(max(abs(D3e)))]

%save coupling matrices
A31_MARSF = A31;
B30_MARSF = B30;
B31_MARSF = B31;

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


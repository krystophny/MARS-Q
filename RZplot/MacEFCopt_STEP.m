% EFC optimization based on various criteria
% Criterion-A: cancellation of 3/1 resonant field (with resistive plasma response)
% Criterion-B: minimize net JxB torque integrated over the plasma volume  
% Criterion-C: SVD based overlap criterion, by performing SVD for the dataset matrix P for deltaB

function MacEFCopt

global SDIRW 

ksave_svd = 0;
ksave_efc = 0;

SDIR  = '/cscratch/liuy/WorkEFC/';
SDIRP = [SDIR 'Database/ResponseMatrix/'];   %folder for storing dataset P = response matrix
SDIRE = [SDIR 'Database/ErrorFields/'];      %folder for storing dataset E = error fields
SDIRC = [SDIR 'Database/CorrectionFields/']; %folder for storing dataset C = correction fields
SDIRD = [SDIR 'Database/DataSave/'];         %folder for storing other needed data files
SDIRR = [SDIR 'Results/'];                   %folder for storing optimization results
SDIRW = [SDIR 'Work/'];                      %working folder

% Step-1: read in response matrix P from database-P
load([SDIRP 'AM_tot.mat']);  %==> AM

% Step-2: read in vacuum EF source in the format of INPUT_BNM
%         from database-E
load([SDIRE 'EM.mat']);  %==> EM

% Step-3: read in vacuum correction field in the format of INPUT_BNM
%         from database-C
load([SDIRC 'CM.mat']); %==> CM

IC = 1.0;  % 1 kAt assumed for correction fields in ERGOS calculations

% cope files to working directory
copyfile([SDIRD 'BPLASMA_SAMPLE'],[SDIRW 'BPLASMA_SAMPLE'],'f');
copyfile([SDIRD 'PROFEQ.OUT'],   [SDIRW 'PROFEQ.OUT'],'f');
copyfile([SDIRD 'RMZM_F_EQAC'],  [SDIRW 'RMZM_F_EQAC'],'f');
copyfile([SDIRD 'RMZM_F_PEST'],  [SDIRW 'RMZM_F_PEST'],'f');

NM = size(EM,1);
NE = size(EM,2);
NC = size(CM,2);
%NE = 1;
%NC = 1;

ICoptA = zeros(NE,NC);
ICoptB = ICoptA;
ICoptC = ICoptA;

BresE31 = zeros(NE,1);
Emax    = BresE31;

if ksave_efc==1
   BresC31 = zeros(NC,1);
   Cmax    = BresC31;
end

% Step-4: perform SVD decomposition for the response matrix, 
%         and save the first (most sensitive) component
if ksave_svd==1
X  = load([SDIRW 'RMZM_F_EQAC']);
NP = round(X(1,2));
X  = load([SDIRW 'BPLASMA_SAMPLE']);
NT = round(X(1,2));

% Get II=row-numbers in P corresponding to B1M inside plasma
II = ones(NP*NM,1);
for m=1:NM
    II((m-1)*NP+[1:NP]) = (m-1)*NT + [1:NP];
end
[U,S,V] = svd(AM(II,:));
Vmax    = V(1,:); Vmax = Vmax(:);
save([SDIRP 'Vmax_tot.mat'],'Vmax');
else
load([SDIRP 'Vmax_tot.mat']);
end

% Step-5.1: compute BresE31 and Emax data 
for j=1:NE
    RE = AM*EM(:,j);
    MacSaveBPLASMA(RE,1);
    BresE31(j) = MacRfaCtBn2;
    Emax(j)    = sum(Vmax.*EM(:,j));
end

% Step-5.2: compute BresC31 and Cmax data 
if ksave_efc==1
for k=1:NC
    RC = AM*CM(:,k);
    MacSaveBPLASMA(RC,1);
    BresC31(k) = MacRfaCtBn2;
    Cmax(k)    = sum(Vmax.*CM(:,k));
end
save([SDIRC 'BresC31_tot.mat'],'BresC31');
save([SDIRC 'Cmax_tot.mat'],'Cmax');
else
load([SDIRC 'BresC31_tot.mat']);
load([SDIRC 'Cmax_tot.mat']);
end

for j=1:NE
disp(['EF number j=',int2str(j),' out of total NE=',int2str(NE)])
for k=1:NC
    % Step-6: EFC optimization with Criterion-A
    ICoptA(j,k) =-IC*BresE31(j)/BresC31(k);

    % Step-7: EFC optimization with Criterion-B
    % Step-7.1: compute plasma response fields
    RE = AM*EM(:,j);
    RC = AM*CM(:,k);
    
    [J1ME,J2ME,B1ME,B2ME] = MacSaveBPLASMA(RE,0);
    [J1MC,J2MC,B1MC,B2MC] = MacSaveBPLASMA(RC,0);

    % Step-7.2: compute mutual net torque coefficients (a,b,c,d)
    %           the net torque with EF + EFC=Tnet(z)= Re[a + b*z + c*z^* + d*|z|^2]
    [a,b,c,d] = MacGetTorqCoeff(J1ME,J2ME,B1ME,B2ME,J1MC,J2MC,B1MC,B2MC);

    % Step-7.3: find analytic solution of z=ICoptB by minimizing Tnet(z)
    ICoptB(j,k)   = -IC*(conj(b)+c)/2/real(d);
    %Tnet_opt = real(a)-0.25/real(d)*(abs(b)^2+2*real(b*c)+abs(c)^2);

    % Step-8: EFC optimization with Criterion-C
    ICoptC(j,k)  = -IC*Emax(j)/Cmax(k);
end
end

% save optimal correction currents
save([SDIRR  'RES_ICoptA.mat'],'ICoptA');
save([SDIRR  'RES_ICoptB.mat'],'ICoptB');
save([SDIRR  'RES_ICoptC.mat'],'ICoptC');
 

%============================================================================
% function for calculating jxb torque coefficients 
function [a,b,c,d] = MacGetTorqCoeff(J1ME,J2ME,B1ME,B2ME,J1MC,J2MC,B1MC,B2MC)

global Mac

N  = Mac.Ns1-1;
II = 1:N;
J2MEN = (J2ME(1:N,:)+J2ME(2:N+1,:))/2;
J2MCN = (J2MC(1:N,:)+J2MC(2:N+1,:))/2;
B1MEN = (B1ME(1:N,:)+B1ME(2:N+1,:))/2;
B1MCN = (B1MC(1:N,:)+B1MC(2:N+1,:))/2;

aa = zeros(size(N,1),1); bb=aa; cc=aa; dd=aa;
for m=1:size(J1ME,2)
    aa = aa + J1ME(II,m).*conj(B2ME(II,m)) - J2MEN(:,m).*conj(B1MEN(:,m));
    bb = bb + J1MC(II,m).*conj(B2ME(II,m)) - J2MCN(:,m).*conj(B1MEN(:,m));
    cc = cc + J1ME(II,m).*conj(B2MC(II,m)) - J2MEN(:,m).*conj(B1MCN(:,m));
    dd = dd + J1MC(II,m).*conj(B2MC(II,m)) - J2MCN(:,m).*conj(B1MCN(:,m));
end

hs = Mac.s(2:N+1)-Mac.s(1:N);
a  = sum(aa.*hs)*2*pi^2;
b  = sum(bb.*hs)*2*pi^2;
c  = sum(cc.*hs)*2*pi^2;
d  = sum(dd.*hs)*2*pi^2;
%============================================================================

%============================================================================
% function for getting and storing JPLASMA and BPLASMA
function [J1M,J2M,B1M,B2M]=MacSaveBPLASMA(RJB,ksave)

global SDIRW

BS = load([SDIRW 'BPLASMA_SAMPLE']);

Mmax = BS(1,1);
Nmax = BS(1,2);
BH   = BS(1:Mmax+1,:);  %header data

RJB = reshape(RJB,Mmax*Nmax,6);
B1M = RJB(:,1);
B2M = RJB(:,2);
B3M = RJB(:,3);
J1M = RJB(:,4);
J2M = RJB(:,5);

if ksave==1
   BN = [BH; real(B1M) imag(B1M) real(B2M) imag(B2M) real(B3M) imag(B3M)];
   save([SDIRW 'BPLASMA.OUT'],'BN','-ascii');
end
%============================================================================








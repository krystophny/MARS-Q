% read in ?PLASMA.OUT files 
% generate (X,Y)-solutions in MARS-F format
% and sane to MARSQ.SOL file, to be read with KSOLTEST=1

SDIR = '/cscratch/liuy/WorkIR_D/BU/';

M0 = 34;

NXC = 6; 
NYC = 8;

KXV1 = 1; 
KXB1 = 2;
KXJ2U= 3; 
KXJ3 = 4; 
KXJ2L= 5; 
KXX1 = 6;

KYV2 = 1;
KYV3 = 6; 
KYB2 = 2; 
KYB3 = 3; 
KYJ1 = 4; 
KYPR = 5; 
KYX2 = 7; 
KYRHOP= 8;

d = load([SDIR 'BPLASMA.OUT']);
MSMAX = round(d(1,1));
ND    = round(d(1,2))-1;
MD    = MSMAX*NXC;
MDY   = MSMAX*NYC;

X = zeros(MD,ND+1);
Y = zeros(MDY,ND+1);

d = d(1+MSMAX+1:end,:);
B1U = d(:,1) + d(:,2)*i;
B2U = d(:,3) + d(:,4)*i;
B3U = d(:,5) + d(:,6)*i;
B1U = reshape(B1U,ND+1,MSMAX);
B2U = reshape(B2U,ND+1,MSMAX);
B3U = reshape(B3U,ND+1,MSMAX);

d = load([SDIR 'VPLASMA.OUT']);
NR = round(d(1,2))-1;
T  = d(1+MSMAX+1:1+MSMAX+NR+1,4);
d = d(1+MSMAX+NR+1+1:end,:);
V1U = d(:,1) + d(:,2)*i;
V2U = d(:,3) + d(:,4)*i;
V3U = d(:,5) + d(:,6)*i;
V1U = reshape(V1U,NR+1,MSMAX);
V2U = reshape(V2U,NR+1,MSMAX);
V3U = reshape(V3U,NR+1,MSMAX);

d = load([SDIR 'XPLASMA_NEW.OUT']);
d = d(1+MSMAX+NR+1+1:end,:);
X1U = d(:,1) + d(:,2)*i;
X2U = d(:,3) + d(:,4)*i;
X3U = d(:,5) + d(:,6)*i;
X1U = reshape(X1U,NR+1,MSMAX);
X2U = reshape(X2U,NR+1,MSMAX);
X3U = reshape(X3U,NR+1,MSMAX);

d = load([SDIR 'JPLASMA_NEW.OUT']);
d = d(1+MSMAX+1:end,:);
J1U = d(:,1) + d(:,2)*i;
J2U = d(:,3) + d(:,4)*i;
J3U = d(:,5) + d(:,6)*i;
J1U = reshape(J1U,ND+1,MSMAX);
J2U = reshape(J2U,ND+1,MSMAX);
J3U = reshape(J3U,ND+1,MSMAX);

d = load([SDIR 'DPLASMA.OUT']);
d = d(1+MSMAX+1:end,:);
RHOP = d(:,1) + d(:,2)*i;
PRE  = d(:,3) + d(:,4)*i;
RHOP = reshape(RHOP,NR+1,MSMAX);
PRE  = reshape(PRE,NR+1,MSMAX);

for MS=1:MSMAX
    LX = (MS-1)*NXC;
    LY = (MS-1)*NYC;

    X(LX + KXV1,1:NR+1) = V1U(:,MS);
    X(LX + KXX1,1:NR+1) = X1U(:,MS);
    X(LX + KXB1,1:ND+1) = B1U(:,MS);
    X(LX + KXJ2U,1:ND+1)= J2U(:,MS);
    X(LX + KXJ3,1:ND+1) = J3U(:,MS);

    if MS==M0
       X(LX + KXB1,1:NR+1) = B1U(1:NR+1,MS)./T;
       X(LX + KXB1,NR+2:ND+1) = B1U(NR+2:ND+1,MS)/T(end);
    end

    Y(LY + KYX2,1:NR)  = X2U(1:NR,MS);
    Y(LY + KYV2,1:NR)  = V2U(1:NR,MS);
    Y(LY + KYV3,1:NR)  = V3U(1:NR,MS);
    Y(LY + KYB2,1:ND)  = B2U(1:ND,MS);
    Y(LY + KYB3,1:ND)  = B3U(1:ND,MS);
    Y(LY + KYJ1,1:ND)  = J1U(1:ND,MS);
    Y(LY + KYPR,1:NR)  = PRE(1:NR,MS);
    Y(LY + KYRHOP,1:NR)= RHOP(1:NR,MS);
end

RES = [0 0 0 1 0 0 0 0 0];
save tmp1 RES -ascii -double

RES = zeros(NR+1,8);
save tmp2 RES -ascii -double

RES = [real(X(:)) imag(X(:)); real(Y(:)) imag(Y(:))];
save tmp3 RES -ascii -double

!cat tmp1 tmp2 tmp3 > MARSQ_NEW.SOL
movefile('MARSQ_NEW.SOL',[SDIR 'MARSQ_NEW.SOL']);

















 



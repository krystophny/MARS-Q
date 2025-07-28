% linear superposition of response data from 3 rows of RMP coils
% for BPLASMA and XPLASMA
% 1. extract b1res and xi_X from individual row response
%    b1res is defined as dimensionless quantity without the q-factor
% 2. 2D scan of coil phasing based upon FOMs from 1
%    specify target FOM values, and scale current amplitude to obtain robust optimium
%    KirkNF15 => b1res_crit=0.5e-4 and xi_X=1.5 mm 
%    ITER 15MA/5.3T case, JOREK==>60kAt for ELM suppression with n=3 ==> xi_X=10mm
% 3. identify best/worst coil phasing from 2, and combine BPLASMA/XPLASMA to plot best/worst response
%    use best/worst coil phasing to also plot corresponding vacuum field
% 4. plot overall optimal results
%
% krot = 0:vac; 1:1e-3; 2:1e-2; 3:5e-2
% kmod = 1:fld; 2:kin

global SDIR n_GLOBAL FIN_RES_b1res FIN_RES_xinX 

SDIR1 = '/home/liuy/Work/STEP-P/Case5/DataN/';
SCOIL = 'UML';
AA    = 1;
R0EXP = [3.6e+3 3.6e+3 3.6e+3 3.6e+3 4.3e+3]; %[mm]
kdb   = 1; %plot threshold database from present devices, valid for krun=4

%SCASE = 'vac';
SCASE = 'rot_5e-2_fld';
n     = 1;
krun  = 4; %see the four steps described above
kbest = 1; %1-best, 0-worst; in action when krun=3

if strcmp(SCASE(1:3),'vac') 
   krot = 0;
   kmod = 1;
else
   if strcmp(SCASE(5:8),'1e-3'), krot=1; end
   if strcmp(SCASE(5:8),'1e-2'), krot=2; end
   if strcmp(SCASE(5:8),'5e-2'), krot=3; end

   if strcmp(SCASE(10:12),'fld'), kmod=1; end
   if strcmp(SCASE(10:12),'kin'), kmod=2; end
end

if krun==1, AA=1;  end
if krun==2, AA=10; end

%result from krun=1
%assume coil current amplitude AA=1kAt
data_b1res=[
%n krot kmod U                          M                          L
1 1    1   -7.4662e-07 - 1.4404e-07i   9.7712e-07 - 1.8647e-07i  -6.2585e-07 + 3.7155e-07i
2 1    1   -2.0760e-07 + 1.7010e-07i   4.6232e-07 - 2.0238e-07i  -2.6434e-07 + 5.0935e-08i
3 1    1   -1.0362e-08 + 1.7112e-09i   3.1441e-08 + 6.6463e-09i  -4.7124e-11 - 8.4930e-09i
4 1    1   -1.0218e-08 + 1.0130e-08i   4.4418e-08 + 2.8982e-09i  -7.6967e-09 - 9.6204e-09i
1 2    1   -7.4125e-07 + 4.5866e-07i   5.4675e-07 - 9.7782e-07i  -1.1382e-07 + 7.0997e-07i
2 2    1    6.0914e-08 + 3.4775e-07i   5.4363e-08 - 6.1574e-07i  -9.3073e-08 + 3.2562e-07i
3 2    1   -5.7148e-09 + 7.6276e-09i   2.3931e-08 - 5.0014e-09i  -8.6157e-09 - 6.0826e-09i
4 2    1   -8.9352e-09 + 8.0668e-09i   4.1086e-08 + 3.4331e-09i  -8.6625e-09 - 9.7917e-09i
1 3    1   -1.6182e-07 + 4.6939e-07i  -2.1189e-07 - 6.2446e-07i   2.3832e-07 + 2.5056e-07i
2 3    1    1.0153e-07 + 2.2208e-07i  -8.1917e-08 - 4.1322e-07i   8.8898e-10 + 2.3647e-07i
3 3    1   -9.4950e-09 + 1.0723e-08i   3.1668e-08 - 4.7138e-10i  -9.2571e-09 - 1.0478e-08i
4 3    1   -7.1975e-09 + 8.8314e-09i   3.8289e-08 + 5.2635e-09i  -6.8728e-09 - 1.1174e-08i
];
data_xin_x=[
%n krot kmod  U                          M                          L
 1 1    1    5.4764e+00 - 8.6995e+00i  -1.4131e+01 + 1.9130e+01i  -3.9975e+00 - 2.1097e+01i
 2 1    1    3.0094e+00 - 8.9558e+00i  -1.0382e+01 + 1.3709e+01i   7.3885e+00 - 5.9802e+00i
 3 1    1   -4.5677e-01 - 6.6896e-01i  -4.4733e-01 + 8.3063e-01i   7.2796e-01 + 9.0987e-02i 
 4 1    1   -1.7338e-01 - 2.6091e-01i  -3.6674e-01 - 9.1488e-03i  -1.0900e-01 + 3.1719e-01i 
 1 2    1    1.1324e+00 + 2.1201e+00i  -3.5961e+00 - 1.1634e+00i   1.6880e+00 - 5.3956e-01i
 2 2    1   -1.2406e+00 - 6.9475e-01i   1.6678e+00 + 1.5333e+00i  -6.3190e-01 - 1.0362e+00i 
 3 2    1   -8.4328e-02 + 8.4944e-03i   1.1829e-01 + 9.2217e-02i   1.6008e-02 - 6.2944e-02i 
 4 2    1   -1.8641e-02 - 2.2596e-02i  -2.3320e-02 + 1.5953e-02i  -5.5086e-03 + 2.1859e-02i 
 1 3    1    2.1006e+00 + 3.3839e-01i  -2.6504e+00 + 1.5071e+00i   7.5125e-01 - 1.3765e+00i
 2 3    1   -4.9684e-01 + 3.3573e-01i   8.9790e-01 - 4.9026e-01i  -4.6452e-01 + 1.5630e-01i
 3 3    1   -8.6782e-02 + 9.5794e-03i   1.1806e-01 + 1.2050e-01i   2.0401e-02 - 7.6477e-02i 
 4 3    1   -6.4746e-03 - 7.5991e-03i  -7.3444e-03 + 1.5464e-02i   1.1857e-03 + 2.3088e-03i 
];

%result from krun=2
%assume coil current amplitude AA=10kAt
%assume PM=0
data_b1max = [
%n krot kmod PUmax        PLmax           b1max     PUmin        PLmin           b1min
 1 1    1   1.5840e+02   2.0160e+02   2.4828e-05   2.5200e+01   3.3120e+02   1.0035e-07 
 2 1    1   1.9440e+02   1.6560e+02   1.0422e-05   3.6000e+01   3.2760e+02   2.6233e-08 
 3 1    1   2.0160e+02   1.0080e+02   5.1129e-07   2.1600e+01   2.8080e+02   1.3144e-07 
 4 1    1   2.3040e+02   1.3320e+02   7.1215e-07   4.6800e+01   3.1320e+02   1.7813e-07 
 1 2    1   1.5120e+02   2.0160e+02   2.7109e-05   2.9160e+02   7.2000e+01   1.3262e-07 
 2 2    1   1.9440e+02   1.6920e+02   1.3098e-05   3.9600e+01   3.2400e+02   1.0118e-07 
 3 2    1   2.2320e+02   1.3320e+02   4.4522e-07   3.9600e+01   3.1320e+02   4.3823e-08 
 4 2    1   2.2680e+02   1.3680e+02   6.6340e-07   4.6800e+01   3.1680e+02   1.6119e-07 
 1 3    1   1.4400e+02   2.0520e+02   1.5016e-05   2.9160e+02   7.2000e+01   2.8908e-08 
 2 3    1   1.9440e+02   1.6920e+02   9.0189e-06   3.4560e+02   1.8000e+01   1.8941e-08 
 3 3    1   2.2680e+02   1.2960e+02   5.9974e-07   4.6800e+01   3.0960e+02   3.4011e-08 
 4 3    1   2.3760e+02   1.2960e+02   6.3159e-07   5.7600e+01   3.0960e+02   1.4141e-07 
];
data_ximax = [
%n krot kmod PUmax        PLmax        ximax        PUmin        PLmin        ximin
 1 1    1   1.8360e+02   2.2680e+02   5.5535e+02   6.8400e+01   2.1600e+01   7.1584e-01 
 2 1    1   1.9800e+02   1.6560e+02   3.6149e+02   3.5280e+02   1.0800e+01   1.4706e+00 
 3 1    1   2.4120e+02   1.1160e+02   2.4868e+01   1.1160e+02   2.3400e+02   1.5078e-01 
 4 1    1   3.0600e+02   7.2000e+01   1.0155e+01   6.8400e+01   3.0600e+02   9.4045e-02 
 1 2    1   1.3680e+02   2.1600e+02   7.9552e+01   2.9520e+02   6.4800e+01   1.6690e-01 
 3 2    1   1.9440e+02   1.6560e+02   4.9007e+01   3.4560e+02   1.8000e+01   1.7319e-01 
 2 2    1   2.2320e+02   1.1520e+02   2.9967e+00   3.9600e+01   2.9880e+02   7.9897e-03 
 4 2    1   2.7360e+02   3.9600e+01   8.0082e-01   1.4040e+02   1.5120e+02   4.1648e-03 
 1 3    1   1.4040e+02   2.1240e+02   6.7445e+01   2.9160e+02   7.2000e+01   3.8642e-01 
 2 3    1   1.8360e+02   1.6920e+02   2.1126e+01   3.4560e+02   1.4400e+01   1.2692e-01 
 3 3    1   2.3040e+02   1.1880e+02   3.3512e+00   5.0400e+01   3.0240e+02   2.3073e-02 
 4 3    1   2.4480e+02   5.4000e+01   2.9696e-01   6.4800e+01   2.3400e+02   4.5446e-02 
];


SDIR2 = [SDIR1 SCASE '/'];
SN    = int2str(n);
SDIR  = SDIR1;
n_GLOBAL = n;

if strcmp(SCASE(1:3),'vac'), SPV='VAC'; end
if strcmp(SCASE(1:3),'rot'), SPV='PLS'; end


b1res = zeros(1,3);  %b1res defoned without q-factor
xin_x = zeros(1,3);  %Xin-X near lower X-point

% 1. extract b1res and xi_X from individual row response
% set CheckErgos=3 and n=n_GLOBAL in MacRfaCtBn2.m
if krun==1
for k=1:length(SCOIL)
    copyfile([SDIR1 'RMZM_F_n' SN '_EQAC'],[SDIR1 'RMZM_F_EQAC'],'f');
    copyfile([SDIR1 'RMZM_F_n' SN '_PEST'],[SDIR1 'RMZM_F_PEST'],'f');
    copyfile([SDIR1 'PROFEQ_n' SN],        [SDIR1 'PROFEQ.OUT'],'f');
    copyfile([SDIR2 'BPLASMA_n' SN '_' SPV '_' SCOIL(k)],[SDIR1 'BPLASMA.OUT'],'f');

    run MacRfaCtBn2.m
    b1res(k) = FIN_RES_b1res;

    if strcmp(SCASE(1:3),'rot')
    copyfile([SDIR1 'RMZM_F_n' SN '_EQAC'],[SDIR1 'RMZM_F.OUT'],'f');
    copyfile([SDIR1 'PROFEQ_n' SN],        [SDIR1 'PROFEQ.OUT'],'f');
    copyfile([SDIR2 'XPLASMA_n' SN '_PLS_' SCOIL(k)],[SDIR1 'XPLASMA.OUT'],'f');

    run MacMainKSTAR_X.m
    xin_x(k) = FIN_RES_xinX;
    end
end
data_b1res=b1res
data_xin_x=xin_x
end

% 2. 2D scan of coil phasing based upon FOMs from 1
% set CheckErgos=3 in MacRfaCtBn2.m
if krun==2
NU    = 101;
NL    = 101;
kB    = find(data_b1res(:,1)==n & data_b1res(:,2)==krot & data_b1res(:,3)==kmod);
kX    = find(data_xin_x(:,1)==n & data_xin_x(:,2)==krot & data_xin_x(:,3)==kmod);

BU = data_b1res(kB,4);
BM = data_b1res(kB,5);
BL = data_b1res(kB,6);

XU = data_xin_x(kX,4);
XM = data_xin_x(kX,5);
XL = data_xin_x(kX,6);

PU = linspace(0,360,NU);
PL = linspace(0,360,NL);

ePU = exp(i*PU*pi/180);
ePL = exp(i*PL*pi/180);

b1res = zeros(NL,NU);
xin_x = zeros(NL,NU);
for ku=1:NU
for kl=1:NL
    b1res(kl,ku) = AA*(BU*ePU(ku) + BM + BL*ePL(kl));
    xin_x(kl,ku) = AA*(XU*ePU(ku) + XM + XL*ePL(kl));
end
end

B = abs(b1res);
X = abs(xin_x);

[Y1,II] = max(B,[],1);
[Y2,JJ] = max(Y1);
IB1 = II(JJ);
JB1 = JJ;
YB1 = Y2;

[Y1,II] = min(B,[],1);
[Y2,JJ] = min(Y1);
IB2 = II(JJ);
JB2 = JJ;
YB2 = Y2;

[Y1,II] = max(X,[],1);
[Y2,JJ] = max(Y1);
IX1 = II(JJ);
JX1 = JJ;
YX1 = Y2;

[Y1,II] = min(X,[],1);
[Y2,JJ] = min(Y1);
IX2 = II(JJ);
JX2 = JJ;
YX2 = Y2;

data_b1max = [PU(JB1) PL(IB1) YB1 PU(JB2) PL(IB2) YB2]
data_ximax = [PU(JX1) PL(IX1) YX1 PU(JX2) PL(IX2) YX2] 

%add q-factor back to b1res when plotting
  d = load([SDIR1 'PROFEQ_n' SN],'-ascii');
q = d(:,2);
m = floor(max(q)*n);
qres = m/n*1e+4;

hf=figure(1);
pcolor(PU,PL,B*qres), colorbar, hold on, shading interp
plot(PU(JB1),PL(IB1),'b+','LineWidth',2,'MarkerSize',12)
plot(PU(JB2),PL(IB2),'wx','LineWidth',2,'MarkerSize',12)
%contour(PU,PL,[1e-4 1e-4],'k-'), hold on,
axis([0 360 0 360]), 
xlabel('\Phi_U-\Phi_M [degree]','FontSize',18,'FontWeight','Bold')
ylabel('\Phi_L-\Phi_M [degree]','FontSize',18,'FontWeight','Bold')
title('|b^1_{res}|x10^4','FontSize',14)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
colormap(hot)

hf=figure(2);
pcolor(PU,PL,X), colorbar, hold on, shading interp
plot(PU(JX1),PL(IX1),'b+','LineWidth',2,'MarkerSize',12)
plot(PU(JX2),PL(IX2),'wx','LineWidth',2,'MarkerSize',12)
%contour(PU,PL,[1e-4 1e-4],'k-'), hold on,
axis([0 360 0 360]), 
xlabel('\Phi_U-\Phi_M [degree]','FontSize',18,'FontWeight','Bold')
ylabel('\Phi_L-\Phi_M [degree]','FontSize',18,'FontWeight','Bold')
title('|\xi_{n,X}| [mm]','FontSize',14)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
colormap(hot)
end

% 3. identify best/worst coil phasing from 2, and combine BPLASMA/XPLASMA to plot best/worst response
% set CheckErgos=0 in MacRfaCtBn2.m
if krun==3
BU = load([SDIR2 'BPLASMA_n' SN '_' SPV '_U']);
BM = load([SDIR2 'BPLASMA_n' SN '_' SPV '_M']);
BL = load([SDIR2 'BPLASMA_n' SN '_' SPV '_L']);

M   = round(BU(1,1));
II  = [0:2]*2;
BUU = (BU(M+2:end,II+1)+BU(M+2:end,II+2)*i);
BMM = (BM(M+2:end,II+1)+BM(M+2:end,II+2)*i);
BLL = (BL(M+2:end,II+1)+BL(M+2:end,II+2)*i);

kB  = find(data_b1max(:,1)==n & data_b1max(:,2)==krot & data_b1max(:,3)==kmod);
if kbest==1
   PUB = data_b1max(kB,4);
   PLB = data_b1max(kB,5);
   %PUB = data_ximax(kB,4);
   %PLB = data_ximax(kB,5);
else
   PUB = data_b1max(kB,7);
   PLB = data_b1max(kB,8);
   %PUB = data_ximax(kB,7);
   %PLB = data_ximax(kB,8);
end
BNN = AA*(BUU*exp(i*PUB*pi/180) + BMM + BLL*exp(i*PLB*pi/180));
BN  = [BU(1:M+1,:); real(BNN(:,1)) imag(BNN(:,1)) ...
                    real(BNN(:,2)) imag(BNN(:,2)) ...
                    real(BNN(:,3)) imag(BNN(:,3))];

save([SDIR1 'BPLASMA.OUT'],'BN','-ascii');

copyfile([SDIR1 'RMZM_F_n' SN '_EQAC'],[SDIR1 'RMZM_F_EQAC'],'f');
copyfile([SDIR1 'RMZM_F_n' SN '_PEST'],[SDIR1 'RMZM_F_PEST'],'f');
copyfile([SDIR1 'PROFEQ_n' SN],        [SDIR1 'PROFEQ.OUT'],'f');

run MacRfaCtBn2.m
b1res = abs(FIN_RES_b1res)

if strcmp(SCASE(1:3),'rot')
XU = load([SDIR2 'XPLASMA_n' SN '_PLS_U']);
XM = load([SDIR2 'XPLASMA_n' SN '_PLS_M']);
XL = load([SDIR2 'XPLASMA_n' SN '_PLS_L']);

M   = round(XU(1,1));
N   = round(XU(1,2));
II  = [0:2]*2;
XUU = (XU(M+N+2:end,II+1)+XU(M+N+2:end,II+2)*i);
XMM = (XM(M+N+2:end,II+1)+XM(M+N+2:end,II+2)*i);
XLL = (XL(M+N+2:end,II+1)+XL(M+N+2:end,II+2)*i);

kX  = find(data_ximax(:,1)==n & data_ximax(:,2)==krot & data_ximax(:,3)==kmod);
if kbest==1
   PUX = data_ximax(kX,4);
   PLX = data_ximax(kX,5);
else
   PUX = data_ximax(kX,7);
   PLX = data_ximax(kX,8);
end
XNN = AA*(XUU*exp(i*PUX*pi/180) + XMM + XLL*exp(i*PLX*pi/180));
XN  = [XU(1:M+N+1,:); real(XNN(:,1)) imag(XNN(:,1)) ...
                      real(XNN(:,2)) imag(XNN(:,2)) ...  
                      real(XNN(:,3)) imag(XNN(:,3))];
save([SDIR1 'XPLASMA.OUT'],'XN','-ascii');

copyfile([SDIR1 'RMZM_F_n' SN '_EQAC'],[SDIR1 'RMZM_F.OUT'],'f');
copyfile([SDIR1 'PROFEQ_n' SN],        [SDIR1 'PROFEQ.OUT'],'f');

run MacMainKSTAR_X.m
xin_x = abs(FIN_RES_xinX)
end

end

if krun==4
%need to add q-factor back to recover CheckErgos=0 
if kdb==1
   run /home/liuy/Work/RMPthreshold/Work/drawPhysThreshold.m
end

run /home/liuy/Work/STEP-P/Work/drawMemo.m
ccol = 'rbkgy';

for kk=1:5
eval(['data_b1max=data_b1max_case' int2str(kk) ';']);
eval(['data_ximax=data_ximax_case' int2str(kk) ';']);

if kdb==0, cc = ccol(kk); end
SS = 'x+s';
SL = [3 1];
for n=[1 2 3 4]
    %get qres-factor at last rational surface
    qa = data0(2*kk,9);
    m  = floor(qa*n);
    qres = m/n;

    II = find(data_b1max(:,1)==n & data_b1max(:,2)~=0);
    db = data_b1max(II,:); 
    II = find(data_ximax(:,1)==n & data_ximax(:,2)~=0);
    dx = data_ximax(II,:); 
    hf1 = figure(10*n+1);
    if kdb==0 
       hf2 = figure(10*n+2);
    elseif kdb==1 
       hf2 = 1; cc=ccol(n); 
    end
    for k=1:size(db,1)
        krot = db(k,2);
        kmod = db(k,3);
        figure(hf1);
        %if krot>1
        %plot(db(k,4),db(k,5),['b' SS(krot)],'LineWidth',SL(kmod),'MarkerSize',9), hold on,
        plot(db(k,4),db(k,5),[cc SS(krot)],'LineWidth',SL(kmod),'MarkerSize',9), hold on,
        %end 

        krot = dx(k,2);
        kmod = dx(k,3);
        %if krot==3
        %plot(dx(k,4),dx(k,5),['g' SS(krot)],'LineWidth',SL(kmod),'MarkerSize',9), hold on,
        plot(dx(k,4),dx(k,5),[cc SS(krot)],'LineWidth',SL(kmod),'MarkerSize',9), hold on,
        %end

        figure(hf2);
        %if krot>1
        if kdb==0
	  loglog(dx(k,6)/R0EXP(kk),db(k,6)*qres,[cc SS(krot)],'LineWidth',SL(kmod),'MarkerSize',9), hold on,
        elseif kdb==1
	  loglog(dx(k,6)/R0EXP(kk),db(k,6)*qres,[cc 'o'],'LineWidth',SL(kmod),'MarkerSize',9,'MarkerFaceColor',cc), hold on,
        end
        %end
    end
    figure(hf1);
    xlabel('\Phi_U-\Phi_M [degree]','FontSize',18,'FontWeight','Bold')
    ylabel('\Phi_L-\Phi_M [degree]','FontSize',18,'FontWeight','Bold')
    ha=get(hf1,'CurrentAxes'); set(ha,'FontSize',16)
    axis([0 360 0 360])
    %if n==4, axis([0 380 -10 360]), end

    figure(hf2);
    xlabel('\xi_n/R_0 ','FontSize',18,'FontWeight','Bold')
    ylabel('b^1_{res}','FontSize',18,'FontWeight','Bold')
    ha=get(hf2,'CurrentAxes'); set(ha,'FontSize',16)
end
end

if kdb==1, figure(hf2), legend('DIII-D','EAST','KSTAR','AUG','STEP'), end
end





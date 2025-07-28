%function MacMain
%% setup for the circular case

global Mac
global SDIR FEEDI KFEED

format short e

Mac.SS = 'b-';
Mac.plot_RZ = 0;
Mac.plot_shape = 1;
Mac.plot_RZM= 0;
Mac.plot_JACOB0 = 0;
Mac.checkRZ = 0;  % compare RZ from Fourier series with RZ from real space
Mac.checkJACOB = 0; % compare jacobian computed from RZ with that computed in CHEASE
Mac.plot_JACOB = 0;

Mac.RunEQ = 0;
Mac.plot_EQ = 0;

Mac.RunJ = 0;
Mac.plot_JM = 0;
Mac.spline_J1 = 0;
Mac.plot_JP = 0;
Mac.plot_JW = 0;
Mac.plot_JW2 = 0;
Mac.plot_Js = 0;
Mac.plot_Jphi = 0;  %plot Jphi along outboard major radius
Mac.plot_J3U = 0;
Mac.plot_J2J3= 0;
Mac.plot_Jpara = 0;
Mac.save_FEEDJ = 0;
Mac.read_FEEDJ = 0;
Mac.plot_FEEDJ = 0;
Mac.plot_shell = 0;
Mac.plot_shdiv = 0;
Mac.DivJ  = 0;
Mac.plot_DivJ = 0;
Mac.ReduceNs = 1;  %reduce radial points when output J on Gaussian points

Mac.RunB = 0;
Mac.RZmap = 0;
Mac.plot_BM = 0;
Mac.plot_Bn = 0;
Mac.spline_B23 = 2;
Mac.plot_B  = 0;
Mac.plot_BR = 0;
Mac.plot_BW = 0;     %plot 2D B-field at wall
Mac.plot_BC = 0;     %plot real-space Bwr,Bwz,Bwphi along physical theta
Mac.plot_Bs = 0;     %plot MARS-F-space B1,B2,B3 along chi at wall
Mac.plot_BWP = 0;    %plot real-space Bwn,Bwt,Bwphi along physical theta 
Mac.plot_BWM = 0;    %plot harmonics of Bwn, Bwt, Bwphi for physical theta
Mac.plot_BWC = 0;    
Mac.plot_BWR = 2;    %plot (BR,BZ) along Mac.rw-surface
Mac.save_BNORM = 0;   %1 <--> BOVACU01; 2<-->BOVACU02; 3<-->BOVACU03; 4<-->BOVACU04; 
Mac.check_BM = 0;    %compare BM from MARS-F and BM from CARIDDI
Mac.plot_BMcheck = 0;
Mac.CNORM = 1.0;
Mac.DivB  = 0;
Mac.plot_DivB = 0;
Mac.plot_JW3  = 0;
Mac.Kratsurf  = 1;
Mac.Iratsurf  = 100;
Mac.FullSol   = 0;
Mac.FootPrint = 0;
Mac.plot_Tjxb = 0;

Mac.RunEF   = 0;  %1: read error field at one surface in vacuum, for generating equivalent surface current
                  %2: read error field at one surface inside plasma, for testing 
Mac.plot_EF = 0;

Mac.RunV    = 0;
Mac.plot_VM = 0;
Mac.plot_V  = 0;
Mac.plot_Vn = 0;

Mac.RunP    = 0;
Mac.plot_PM = 0;
Mac.plot_PMM= 0;

Mac.plot_coil = 1;
Mac.NGauss = 2;      % order of Gauss quadrature integration
Mac.Nchi   = 513;    % number of points along poloidal angle 'chi'
Mac.Nm2    = 100;     % number of poloidal harmonics for RZ-construction
Mac.Nm0    = [];     % number of poloidal harmonics for equilibrium
Mac.Nm1    = [];     % number of poloidal harmonics for stability
Mac.Ns1    = [];     % number of radial points in plasmas
Mac.Ns2    = [];     % number of radial points in vacuum
Mac.Ns     = [];     % = Mac.Ns1 + MacNs2;
Mac.chi    = [];     % poloidal coordinate 'chi'
Mac.phi    = [];     % toroidal coordinate 'phi'
Mac.s      = [];     % radial coordinate 's'
Mac.rw     = [];     % radial position of flux surfaces to be saved
Mac.Mm     = [];     % array of poloidal harmonics for stability
Mac.n      = [];     % toroidal mode number
Mac.Norm   = 1;  % Apply normalization
Mac.R0EXP  = 2.0;    % active if Mac.Norm=0 
Mac.B0EXP  = 1.0;    % active if Mac.Norm=0 
Mac.HWEXP  = 0.01;   % wall thickness [m]
Mac.coil   = [1.24 0.1334 0.2584; 1.24 -0.2530 -0.123];     
           % coils position defined as [rc chi_1 chi_2]
Mac.resetCoil = 1;
Mac.coilN  = [2.0084e+00 3.9247e-01 2.1527e+00 1.7257e-01;
                    2.1527e+00  -1.7257e-01 2.0084e+00 -3.9247e-01]; 
                 % HL-2A ELM coils: [R1 Z1 R2 Z2;...]
Mac.rw     = [1.0 1.368];
Mac.rs     = Mac.rw(2); % the position of equivalent sheet current from CARIDDI
Mac.cut_RZ = -1;
Mac.plot_thickness = 0;
Mac.CHIS   = 0;


phas = 0.0*pi;
Mac.ALPHA   = 1.0;                     
Mac.ShellCurrNorm = 1.0;
Mac.BNORM = 1.0*exp(phas*i);
Mac.JNORM = 1.0*exp(phas*i);
Mac.PNORM = 1.0*exp(phas*i); 
Mac.VNORM = 1.0*exp(phas*i);

SDIR = '/home/yliu/HL-2A/SHOT29676/Work/';

if Mac.checkRZ, Mac.plot=1; end
Mac.chi = linspace(-pi,pi,Mac.Nchi);
Mac.phi = linspace(0,2*pi,101)';

% RZ coordinates
[RM,ZM] = MacReadRMZM([SDIR 'RMZM_F_n1_EQAC']);
[R,Z] = MacGetRZ(RM,ZM);
if Mac.checkRZ, MacCheckRZ('../TestCirc/RMZM_R'), end

% get R(s) along mid-plane inside the plasma
% better to use GEOM coordinates for this
if 1==0
II   = round((Mac.Nchi+1)/2); III = 1:Mac.Ns1;
res = [Mac.s(III) R(III,1)*Mac.R0EXP R(III,II)*Mac.R0EXP];
save R_VS_S_MIDPLANE res -ascii -double
plot(res(:,2),Z(III,1)*Mac.R0EXP,'b-','LineWidth',3), hold on,
plot(res(:,3),Z(III,II)*Mac.R0EXP,'b-','LineWidth',3), hold on,
end

[Y,Imid] = max(R(Mac.Ns1,:))
[Y,Ixpt] = min(Z(Mac.Ns1,:))

if Mac.RunEQ
[QEQ] = MacReadPROFEQ([SDIR 'PROFEQ.OUT']);
end

% normalization for B-field
fac_c = 0.3852;      %=I_n/I_0
I0EXP = 1.0e+3*fac_c;      %[kA] vacuum shot
mu0   = 4e-7*pi;
fac   = mu0/Mac.R0EXP;
FEEDI = 1.0;
%Mac.BNORM   = Mac.B0EXP;
%Mac.BNORM  = Mac.BNORM*fac*I0EXP/FEEDI*1e+4;  %in Gauss/kA
%Mac.VNORM  = Mac.BNORM*Mac.R0EXP/Mac.B0EXP/1e+4; %in meter/kA
%Mac.JNORM  = Mac.BNORM*(Mac.B0EXP/Mac.R0EXP/4.0e-7/pi)/Mac.B0EXP/1e+4/1e+3; %in kA/kA
if Mac.RunEF>0
   %Mac.BNORM   = Mac.B0EXP*1e+4*(4.8881e-04+1.4664e-03i)/1.05;
   %Mac.BNORM   = Mac.B0EXP*1e+4*(9.7775e-04+9.7775e-04i);
      Mac.BNORM   = Mac.B0EXP*1e+4;
   %Mac.BNORM   = 1.0;
end

% RZ coordinates for the coils
if Mac.resetCoil, MacResetCoil(R,Z), end
if length(Mac.coil)>0, [Rc,Zc] = MacGetRZcoil(RM,ZM); end

% unit vectors e_s and e_chi and jacobian
[dRds,dZds,dRdchi,dZdchi,jacobian] = MacGetUnitVec(R,Z);
if Mac.checkJACOB, MacCheckJacob([SDIR 'JACOBIAN_fix'],jacobian), end

%read in error field data
if Mac.RunEF==1
   MacGetB1EF([SDIR 'surfs_ef_2007554.dat'],R,Z,dRdchi,dZdchi);
elseif Mac.RunEF==2
   %MacGetB1EF([SDIR 'surft_ef_2007554.dat'],R,Z,dRdchi,dZdchi);
   MacGetB1EF([SDIR 'surfs_ef_2007554.dat'],R,Z,dRdchi,dZdchi);
end

% perturbed plasma current
if Mac.RunJ
[JM1,JM2,JM3] = MacReadJPLASMA([SDIR 'JPLASMA.OUT']);
if Mac.DivJ & Mac.spline_J1==0, MacCheckDivJ(JM1, JM2, JM3); end
if Mac.plot_J3U, MacGetJ3UMAX('../SHOT117327/J3UMAX_RFA_I2P05E0T2e4W8em7'), end
if 1==1
[J1,J2,J3] = MacGetJ123(JM1,JM2,JM3,R,Z);
[Jr, Jz, Jphi] = MacGetJphysC(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,J1,J2,J3);
MacGetJatGaussQuad(R, Z, Jr, Jz, Jphi);
MacGetJatWall(R, Z, Jr, Jz, Jphi);
[Jrho, Jchi, Jphi] = MacGetJphysT(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,J1,J2,J3);
end

% equivalent sheet current from CARIDDI
if length(Mac.rs)==1 & Mac.read_FEEDJ 
[TC, JrC, JzC, JphiC] = MacReadJCARIDDI('../TestCirc/shell_F_iteration_19.asc');
[CHIG,WG,J1G,J2G,J3G] = MacGetJ123G(TC,JrC,JzC,JphiC,R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
[JM1C,JM2C,JM3C] = MacGetJM123(CHIG,WG,J1G,J2G,J3G);
end

end

% perturbed magnetic field
if Mac.RunB
%MacReadGAMMA('../TestCirc/GAMMA_TMP');
[BM1,BM2,BM3] = MacReadBPLASMA([SDIR 'BPLASMA_PLS']);

if 1==1

[B1,B2,B3,Bn] = MacGetB123(BM1,BM2,BM3,R,dRdchi,dZdchi);
%MacGetBsn(B1,B2);
if Mac.DivB & Mac.spline_B23==0, MacCheckDivB(BM1, BM2, BM3); end
%[Brho, Bchi, Bphi] = MacGetBphysT(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,B1,B2,B3);
[Br, Bz, Bphi] = MacGetBphysC(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,B1,B2,B3);
%[Bnm,Btm,Bpm]=MacGetBatWall(R, Z, Br, Bz, Bphi);
%MacGetBnWall(R, Z, Br, Bz, Bphi);

if Mac.FootPrint>0, MacGetBforFootPrint(R,Z,Br,Bz); end

if length(Mac.rw) == 1 & Mac.save_BNORM == 1
[TC, BrC, BzC, BphiC] = MacReadBCARIDDI('../TestCirc/MacDataBs_CARIDDI');
[TC, BrC, BzC, BphiC] = MacReadBMARSF('../TestCirc/MacDataBs_MARSF');
[CHIG,WG,B1G,B2G,B3G] = MacGetB123G(TC,BrC,BzC,BphiC,R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
[BM1C,BM2C,BM3C] = MacGetBM123(CHIG,WG,B1G,B2G,B3G);
if Mac.check_BM, MacCheckBM(BM1,BM2,BM3,-BM1C,-BM2C,-BM3C); end
elseif length(Mac.rw) == 2 & (Mac.save_BNORM == 3|Mac.save_BNORM == 4)
MacSplitBCARIDDI('../TestCirc/MacDataBs_2N'); 
rw_save = Mac.rw;
Mac.rw = rw_save(1);
[TC, BrC, BzC, BphiC] = MacReadBCARIDDI('MacDataBs_TMP1');
[CHIG,WG,B1G,B2G,B3G] = MacGetB123G(TC,BrC,BzC,BphiC,R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
[BM1C1,BM2C1,BM3C1] = MacGetBM123(CHIG,WG,B1G,B2G,B3G);
Mac.rw = rw_save(2);
[TC, BrC, BzC, BphiC] = MacReadBCARIDDI('MacDataBs_TMP2');
[CHIG,WG,B1G,B2G,B3G] = MacGetB123G(TC,BrC,BzC,BphiC,R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
[BM1C2,BM2C2,BM3C2] = MacGetBM123(CHIG,WG,B1G,B2G,B3G);
Mac.rw = rw_save;
tmp = BM1C1;
tmp1 = [real(tmp(:)) imag(tmp(:))];
tmp = BM2C2;
tmp2 = [real(tmp(:)) imag(tmp(:))];
if Mac.save_BNORM == 3
  tmp = BM1C1./BM2C2;
elseif Mac.save_BNORM == 4
  tmp = BM1C1/sum(BM2C2);
end
tmp3 = [real(tmp(:)) imag(tmp(:))];
tmp = [tmp1 tmp2 tmp3];
if Mac.save_BNORM == 3
  save BNORM03_CARIDDI tmp3 -ascii -double
elseif Mac.save_BNORM == 4
  save BNORM04_CARIDDI tmp3 -ascii -double
end
end
end
end

% compute toroidal component of the jxb torque density
if Mac.RunJ & Mac.RunB & Mac.spline_B23==0 & Mac.spline_J1==0 
[Tjxb,Ajxb] = MacGetTorqJXB(BM1,BM2,JM1,JM2,jacobian);
res_JXB = Ajxb(end)*Mac.R0EXP^3*Mac.B0EXP^2/4e-7/pi 
end

% perturbed velocity
if Mac.RunV
[VM1,VM2,VM3,DPSIDS,T] = MacReadVPLASMA([SDIR 'XPLASMA.OUT']);
[V1,V2,V3,Vn] = MacGetV123(VM1,VM2,VM3,R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,DPSIDS,T);
%[Vr, Vz, Vphi] = MacGetVphys(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,V1,V2,V3);
%data = load([SDIR 'PROFEQ.OUT']); density = data(:,5);
%den_pump_out   = MacGetDenPumpOut(R,V1,jacobian,density)

[Y,Imid] = max(R(Mac.Ns1,:))
[Y,Ixpt] = min(Z(Mac.Ns1,:))
res_Vn = [max(abs(Vn(end,:))) mean(abs(Vn(end,:))) abs(Vn(end,Imid)) abs(Vn(end,Ixpt))]*Mac.R0EXP*1e+3 
end


% perturbed pressure
if Mac.RunP
[PM1,PM2,PM3] = MacReadPPLASMA([SDIR 'PPLASMA.OUT']);
end

if 1==0
%compare B1 and V1
data = load([SDIR 'PROFEQ_PEST']); 
psip = data(end,12)
F    = data(end,13)
N    = Mac.Ns1;
B1E  = B1(N,:);
V1E  = V1(N,:);
RE   = R(N,:);
JE   = jacobian(N,:);
ntor = -3;
x    = Mac.chi;
xx   = (x(1:end-1)+x(2:end))/2;
yy   = diff(V1E)./diff(x);
zz   = spline(xx,yy,x);

Y1E  = psip*zz + JE*F./RE.^2*i*ntor.*V1E;

figure
plot(x,JE)
figure
plot(x,RE)
figure
plot(x,abs(V1E))
figure
plot(x,abs(zz))

figure
plot(x,real(B1E),'r-','LineWidth',2), hold on,
plot(x,real(Y1E),'b--','LineWidth',2), hold on,

figure
plot(x,imag(B1E),'r-','LineWidth',2), hold on,
plot(x,imag(Y1E),'b--','LineWidth',2), hold on,
end

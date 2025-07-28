%function MacMain

global Mac
global SDIR FEEDI KFEED

format short e

Mac.SS = 'b-';
Mac.plot_RZ     = 0;
Mac.plot_shape  = 1;
Mac.plot_thickness = 0;
Mac.plot_RZM    = 0;
Mac.checkRZ     = 0;  % compare RZ from Fourier series with RZ from real space
Mac.checkJACOB  = 0; % compare jacobian computed from RZ with that computed in CHEASE
Mac.plot_JACOB  = 0;
Mac.plot_JACOB0 = 0;

Mac.RunEQ = 0;
Mac.plot_EQ1 = 0;
Mac.plot_EQ2 = 0;
Mac.plot_EQ3 = 0;

Mac.RunJ = 0;
Mac.plot_JM = 1;
Mac.spline_J1 = 1;
Mac.plot_JP = 0;
Mac.plot_JW = 0;
Mac.plot_JW2 = 0;
Mac.plot_Js = 0;
Mac.plot_J3U = 0;
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
Mac.spline_B23 = 2;
Mac.plot_B  = 0;
Mac.plot_BR = 0;
Mac.plot_BW = 0;     %plot 2D B-field at wall
Mac.plot_BC = 0;     %plot real-space Bwr,Bwz,Bwphi along physical theta
Mac.plot_Bs = 0;     %plot MARS-F-space B1,B2,B3 along chi at wall
Mac.plot_BWP = 0;    %plot real-space Bwn,Bwt,Bwphi along physical theta 
Mac.plot_BWC = 0;    %2D plot of Re(Bwn) along theta and phi
Mac.plot_BWM = 0;    %plot harmonics of Bwn, Bwt, Bwphi for physical theta
Mac.save_BNORM = 0;   %1 <--> BOVACU01; 2<-->BOVACU02; 3<-->BOVACU03; 4<-->BOVACU04; 
Mac.check_BM = 0;    %compare BM from MARS-F and BM from CARIDDI
Mac.plot_BMcheck = 0;
Mac.CNORM = 1.0;
Mac.DivB  = 0;
Mac.plot_DivB = 0;

%plasma velocity perturbation
Mac.RunV = 0;
Mac.plot_VM = 3;
Mac.plot_V = 0;

%compute the NTV torque
Mac.RunNTV = 0;        %compute NTV torque due to B-field perturbation
Mac.dataEXP = 4;      %read dataDEN,dataTEM,dataROT 
Mac.plot_nui    = 0;  %ion collision frequency
Mac.plot_Vti    = 0;  %ion thermal velocity
Mac.plot_BHAM   = 0;  %plot Fourier harmonics of B in Hamada coordinates
Mac.plot_ChiHAM = 0;  %plot Chi(theta) for Hamada coodinates at plasma boundary 
Mac.plot_NTVR   = 0;  %plot Rmid & Rave for NTV calculations
Mac.plot_NTVBeq = 0;  %plot Rmid & Rave for NTV calculations
Mac.plot_NTVW = 0;    %plot Wnm for NTV calculations
Mac.plot_NTVCnu = 0; %plot ion collisionality Cnu
Mac.plot_NTVSv = 0;  %plot volume of innertial layer
Mac.plot_NTV  = 17;   %plot the NTV torq

%compute the jxb torque
Mac.RunJXB = 0;
Mac.plot_JXB = 0;

%plot total torque
Mac.plot_Torq = 0;

%basic parameter setting
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
Mac.Norm   = 1;      % Apply normalization
Mac.R0EXP  = 2.0;    % active if Mac.Norm=0 
Mac.B0EXP  = 1.0;    % active if Mac.Norm=0 
Mac.HWEXP  = 0.01;   % wall thickness [m]
Mac.coil   = [2.9 -0.2596 0.2180; 1.45 -0.1340 0.042];     % coils position defined as [rc chi_1 chi_2]
Mac.resetCoil = 1;
%Mac.coilN  = [3.23 0 1.6];     % C-coil: [Rc Zc W]
%Mac.coilN  = [5.5   3.09 5.5  -2.99;  2.04 -1.05 3.56 -1.09; 4.41 0.64 4.41 -0.64];     
Mac.coilN  = [5.5   3.09 5.5  -2.99;  4.41 0.64 4.41 -0.64];     
%Mac.coilN  = [5.5   3.09 5.5  -2.99;  4.41 0.79 4.41 -0.49];   %sensor shifted upwards by 15cm  
Mac.rw     = [1.0 1.4350];
Mac.rs     = 1.05; % the position of equivalent sheet current from CARIDDI
Mac.cut_RZ = -1;
Mac.CHIS   = 0;
%Mac.cut_RZ      = 1.426;

Mac.ALPHA   = 1.0;                     
Mac.ShellCurrNorm = 1.0;
BPHASE = 0.0;
Mac.BNORM = exp(i*pi*BPHASE);
Mac.JNORM = exp(i*pi*BPHASE);
Mac.VNORM = exp(i*pi*BPHASE);

SDIR = '/.automount/funsrv1/root/home/yliu/JET77329/';

if Mac.checkRZ, Mac.plot=1; end
Mac.chi = linspace(-0.3*pi,0.3*pi,Mac.Nchi);
Mac.phi = linspace(0,2*pi,Mac.Nchi)';

% RZ coordinates
[RM,ZM] = MacReadRMZM([SDIR 'RMZM_F']);
[R,Z] = MacGetRZ(RM,ZM);
if Mac.checkRZ, MacCheckRZ('../TestCirc/RMZM_R'), end

if Mac.RunEQ | Mac.RunNTV
[QEQ] = MacReadPROFEQ([SDIR 'PROFEQ_EQARC']);
end

% RZ coordinates for the coils
if Mac.plot_coil > 0
if Mac.resetCoil, MacResetCoil(R,Z), end
if length(Mac.coil)>0, [Rc,Zc] = MacGetRZcoil(RM,ZM); end
end

% unit vectors e_s and e_chi and jacobian
[dRds,dZds,dRdchi,dZdchi,jacobian] = MacGetUnitVec(R,Z);
if Mac.checkJACOB, MacCheckJacob([SDIR 'JACOB_F_EQARC'],jacobian), end

% normalization for B-field and J-field
%FEEDI = 0.3214;         %BPLASMA_EQARC
%FEEDI = 0.3316;         %BPLASMA_VAC_EQARC
%I0EXP = 0.9e+3*16;      %[A] vacuum shot
%mu0   = 4e-7*pi;
%Mac.BNORM = mu0*I0EXP/Mac.R0EXP/FEEDI;
%Mac.JNORM = I0EXP/FEEDI/Mac.R0EXP^2; 

% perturbed plasma current
if Mac.RunJ | Mac.RunJXB
[JM1,JM2,JM3] = MacReadJPLASMA([SDIR 'JPLASMA_EQARC']);
if Mac.DivJ & Mac.spline_J1==0, MacCheckDivJ(JM1, JM2, JM3); end
if Mac.plot_J3U, MacGetJ3UMAX('../SHOT117327/J3UMAX_RFA_I2P05E0T2e4W8em7'), end
[J1,J2,J3] = MacGetJ123(JM1,JM2,JM3);
[Jr, Jz, Jphi] = MacGetJphysC(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,J1,J2,J3);
MacGetJatGaussQuad(R, Z, Jr, Jz, Jphi);
MacGetJatWall(R, Z, Jr, Jz, Jphi);
[Jrho, Jchi, Jphi] = MacGetJphysT(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,J1,J2,J3);

% equivalent sheet current from CARIDDI
if length(Mac.rs)==1 & Mac.read_FEEDJ 
[TC, JrC, JzC, JphiC] = MacReadJCARIDDI('../TestCirc/shell_F_iteration_19.asc');
[CHIG,WG,J1G,J2G,J3G] = MacGetJ123G(TC,JrC,JzC,JphiC,R,Z,dRds,dZds,dRdchi,dZdchi,jacobian);
[JM1C,JM2C,JM3C] = MacGetJM123(CHIG,WG,J1G,J2G,J3G);
end

end

% perturbed magnetic field
if Mac.RunB | Mac.RunNTV | Mac.RunJXB
[BM1,BM2,BM3] = MacReadBPLASMA([SDIR 'BPLASMA']);
[B1,B2,B3,Bn] = MacGetB123(BM1,BM2,BM3,R,dRdchi,dZdchi);
%MacGetBsn(B1,B2);
if Mac.DivB & Mac.spline_B23==0, MacCheckDivB(BM1, BM2, BM3); end
%[Brho, Bchi, Bphi] = MacGetBphysT(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,B1,B2,B3);
[Br, Bz, Bphi] = MacGetBphysC(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,B1,B2,B3);
[Bnm,Btm,Bpm]=MacGetBatWall(R, Z, Br, Bz, Bphi);

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


% perturbed velocity
if Mac.RunV
[VM1,VM2,VM3,DPSIDS,T] = MacReadVPLASMA([SDIR 'VPLASMA']);
[V1,V2,V3] = MacGetV123(VM1,VM2,VM3,R,dRds,dZds,dRdchi,dZdchi,jacobian,DPSIDS,T);
[Vr, Vz, Vphi] = MacGetVphys(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,V1,V2,V3);
end


[Rmid,Rave] = MacGetNTVR(R,Z,jacobian);

% compute NTV torque
if Mac.RunNTV
[RMH,ZMH] = MacReadRMZM([SDIR 'RMZM_F_HAMADA']);
[RH,ZH] = MacGetRZ(RMH,ZMH);

BMH = MacGetNTVB(R,Z,RH,ZH,sqrt(abs(Br).^2+abs(Bz).^2+abs(Bphi).^2));
[Bmid,Btmid,BBtave,Btave] = MacGetNTVBeq(R,Z,dRds,dZds,jacobian);
Wnm = MacGetNTVW;
[Vti,nui,rot,peqi] = MacGetNTVV(Rmid);
Sv  = MacGetNTVSv(Rmid,jacobian);
TorqNTV = MacGetNTVTorq(BMH,Rmid,Rave,Vti,Bmid,Btmid,BBtave,Btave,Wnm,nui,rot,peqi,Sv);
end

% compute jxb torque
if Mac.RunJXB
TorqJXB = MacGetJXBTorq(Rmid,JM1,JM2,BM1,BM2);
end

if Mac.plot_Torq>0 & Mac.RunNTV & Mac.RunJXB
Torq = MacGetTotTorq(Rmid,TorqNTV,TorqJXB);
end

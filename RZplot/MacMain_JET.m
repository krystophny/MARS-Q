function MacMain

global Mac

Mac.plot   = false;  % plot figures ?
Mac.plot_B = true;
Mac.checkRZ = false;  % compare RZ from Fourier series with RZ from real space
Mac.checkJACOB = false; % compare jacobian computed from RZ with that computed in CHEASE
Mac.NGauss = 2;      % order of Gauss quadrature integration
Mac.Nchi   = 101;    % number of points along poloidal angle 'chi'
Mac.Nm2    = 18;     % number of poloidal harmonics for RZ-construction
Mac.Nm0    = [];     % number of poloidal harmonics for equilibrium
Mac.Nm1    = [];     % number of poloidal harmonics for stability
Mac.Ns1    = [];     % number of radial points in plasmas
Mac.Ns2    = [];     % number of radial points in vacuum
Mac.Ns     = [];     % = Mac.Ns1 + MacNs2;
Mac.chi    = [];     % poloidal coordinate 'chi'
Mac.s      = [];     % radial coordinate 's'
Mac.rw     = [];     % radial position of flux surfaces to be saved
Mac.Mm     = [];     % array of poloidal harmonics for stability
Mac.n      = [];     % toroidal mode number
Mac.R0EXP  = [];     % normalization factor for CHEASE and MARS-F

if Mac.checkRZ, Mac.plot=true; end
Mac.chi = linspace(0,2*pi,Mac.Nchi);
Mac.rw  = [1.0 1.325 3.5];

% RZ coordinates
[RM,ZM] = MacReadRMZM('RMZM_FJET');
[R,Z] = MacGetRZ(RM,ZM);
if Mac.checkRZ, MacCheckRZ('RMZM_RJET'), end

% perturbed plasma current
[JM1,JM2,JM3] = MacReadJPLASMA('JPLASMA_JET');
[J1,J2,J3] = MacGetJ123(JM1,JM2,JM3);
[Jr, Jz, Jphi, jacobian] = MacGetJrzphi(R,Z,J1,J2,J3);
MacGetJatGaussQuad(R, Z, Jr, Jz, Jphi)

% perturbed magnetic field
[BM1p,BM2p,BM3p] = MacReadBPLASMA('BPLASMA_JETp');
[BM1m,BM2m,BM3m] = MacReadBPLASMA('BPLASMA_JETm');
BM1 = 0.5*(BM1p+BM1m); BM2 = 0.5*(BM2p+BM2m); BM3 = 0.5*(BM3p+BM3m);
[B1,B2,B3] = MacGetB123(BM1,BM2,BM3);
[Brho, Bchi, Bphi] = MacGetBphys(R,Z,B1,B2,B3);

% check jacobian
if Mac.checkJACOB, MacCheckJacob('JACOB_FJET',jacobian), end


% This memo-file records key information during the samples runs of MARS-*
% based on an artificial but full toroidal equilibrium with aspect ratio of 3.
% YQ Liu, 2014-05-07

% The following sub-folders are created:
% ---------------------------------------------------------
% Memo:       contains this memo-file
% Data:       store key input data for CHEASE and MARS-* runs
% Equi:       computing equilibrium with CHEASE
% KinK:       computing stability of ideal kink mode
% RwmFluid:   computing stability of fluid RWM
% RwmKinetic: computing RWM stability with kinetic effects
% RmpLin:     computing plasma response to RMP fields
% RmpQLin:    modeling RMP field penetration
% Tm:         computing growth rate of TM
% -----------------------------------------------------------

% Equi: computing equilibrium with CHEASE
% input file: EXPEQ_EXAMPLE
% run: Rchease_EXAMPLE_base (for all cases below except Tm)
%      fix total plasma current by setting NCSCAL=2,CURRT=0.30
%      assume the base case with CFBAL=1.0
% output file: log_chease_base
% record basic plasma parameters
dataEqui = [
%CFBAL R0EXP[m] B0EXP[T] I0EXP[MA] BETAN(GEXP) q0     q95    qa      NW 
 1.0   3.0      1.5      1.0743    1.2262      1.4792 5.0237 6.2477  23
];
    
% Kink: perform pressure scan for computing the no-wall and ideal-wall
% betan limits, assuming a wall position at rw=1.2308a
% input file: PROFDEN_EXAMPLE.IN
%             OUTRMAR, OUTVMAR
% run: Rmars_EXAMPLE_Kink
% output file: RESULT.OUT, log_mars
% record pressure scan results for the n=1 ideal kink growth rate 
% first scan ideal-wall betan limit: set NV=NW-1 in MARS-* runs
%                                    set Im(TALPHA1)=0.0
dataKink_idwall = [
%CFBAL BETAN  q0     qa     NW  Re(TALPHA1) Re(gamma)*tauA  
 5.0   6.2274 1.1899 9.7875 22  5.00e-2     4.31991E-02
 4.8   5.9868 1.2028 9.4365 24  3.50e-2     3.37854E-02
 4.6   5.7446 1.2157 9.1219 31  2.50e-2     2.33546E-02
 4.4   5.5008 1.2285 8.8378 21  1.50e-2     2.05345E-02
 4.3   5.3784 1.2350 8.7059 21  1.60e-2     1.70907E-02
 4.2   5.2556 1.2414 8.5799 22  1.50e-2     1.26812E-02
 4.1   5.1324 1.2479 8.4595 23  0.70e-2     7.98842E-03
 4.0   5.0089 1.2544 8.3443 24  0.60e-2     3.03150E-03
];
data = dataKink_idwall;
x    = data(:,2);
y    = data(:,end);
xx   = linspace(x(1),x(end)-0.2,101);
yy   = spline(x,y,xx);
figure(1)
plot(x,y,'ro',xx,yy,'b-'), hold on,
% ==> betan_nw=4.9369

% now scan no-wall betan limit: set NV=160 in MARS-* runs
%                               set Im(TALPHA1)=0.0
dataKink_nowall = [
%CFBAL BETAN  q0     qa     NW  Re(TALPHA1) Re(gamma)*tauA  
 4.0   5.0089 1.2544 8.3443 24  1.00e-1     1.03939E-01 
 3.5   4.3864 1.2876 7.8343 20  6.00e-2     6.04377E-02
 3.2   4.0096 1.3081 7.5719 21  3.00e-2     3.56569E-02
 3.0   3.7573 1.3220 7.4116 22  1.50e-2     1.89481E-02 
 2.9   3.6308 1.3291 7.3353 23  1.00e-2     1.19395E-02
 2.8   3.5042 1.3362 7.2615 24  0.60e-2     5.64874E-03
];
data = dataKink_nowall;
x    = data(:,2);
y    = data(:,end);
xx   = linspace(x(1),x(end)-0.2,101);
yy   = spline(x,y,xx);
figure(1)
plot(x,y,'rs',xx,yy,'b--'), hold on,
% ==> betan_iw=3.3887

% work out correspodence between CFBAL,Cbeta,BETAN
% between no-wall and ideal-wall limits
beta_nw = 3.3887;
beta_iw = 4.9369;
Cbeta   = linspace(0,1,11)';
betan   = beta_nw + Cbeta*(beta_iw-beta_nw);
data    = [dataKink_idwall(:,1:2); dataKink_nowall(2:end,1:2)];
x       = [data(:,1); 1.0];
y       = [data(:,2); 1.2262];
CFBAL   = spline(y,x,betan);
figure(2)
plot(x,y,'b-',CFBAL,betan,'ro')
% ==> 
res_Kink = [
%Cbeta  CFBAL     BETAN
 0.0    2.7088    3.3887
 0.1    2.8311    3.5435
 0.2    2.9534    3.6983
 0.3    3.0759    3.8532
 0.4    3.1987    4.0080
 0.5    3.3218    4.1628
 0.6    3.4451    4.3176
 0.7    3.5687    4.4724
 0.8    3.6927    4.6273
 0.9    3.8171    4.7821
 1.0    3.9418    4.9369
];
 
% RwmFluid: compute stability of RWM in fluid approximation
% consider as an example Cbeta=0.5 case
% assume a resistive wall at rw=1.2308a ==> NW=21=IWALL
% assume TAUW=1.0e+4
% assume a weak parallel sound wave damping for the RWM: PVISC=0.1
% input file: PROFDEN_EXAMPLE.IN, PROFROT_EXAMPLE.IN
%             OUTRMAR, OUTVMAR
% CASE#1: run <Rmars_EXAMPLE_RwmFluid_walltime> to compute wall time
%         output: RESULT_walltime, log_mars_walltime
%         ==> computed wall time TAUWC=-1/(-8.86148E-05)=1.1285e+04tauA
% CASE#2: run <Rmars_EXAMPLE_RwmFluid_static> to compute RWM growth rate w/o flow
%         output: RESULT_static, log_mars_static
%         ==> computed growth rate = gamma*TAUWC=7.2903
% CASE#3: run <Rmars_EXAMPLE_RwmFluid_flow> to compute RWM growth rate w/ flow
%         output: RESULT_flow, log_mars_flow
%         ==> computed growth rate = gamma*TAUWC=7.2817+0.1709i
dataRwmFluid = [
%CASE# CFBAL  BETAN  q0     qa     NW  TALPHA1(Re:Im) gamma*tauA(Re:Im)  
 1     3.3218 4.1628 1.2997 7.6750 21 -2.0e-5 0.0    -8.86148E-05  1.25206E-09
 2     3.3218 4.1628 1.2997 7.6750 21  5.0e-4 0.0     6.46017E-04  1.71190E-14
 3     3.3218 4.1628 1.2997 7.6750 21  6.5e-4 2.0e-5  6.45251E-04  1.51480E-05
];

% RwmKinetic: compute stability of RWM with drift kinetic effects
% consider an example of non-perturbative run
% consider the same equilibrium as in RwmFluid case
% assume OMEGACI0=40.0
% input file: PROFDEN_EXAMPLE.IN, PROFROT_EXAMPLE.IN, PROFWE_EXAMPLE.IN
%             OUTRMAR, OUTVMAR
%             note that here we simply assume PROFROT=PROFWE
% CASE#1: run <Rmars_EXAMPLE_RwmKinetic_1>: consider only precessional
%         drift resonances of thermal ions and electrons. No energetic
%         particle contribution is included. 
%         output: RESULT_1, ENERGY_1, log_mars_1
%         note: the result here is just one example. There might be more than
%               one unstable branch with the non-perturbative approach. The proper
%               way of performing non-perturbative hybrid runs is to scan 
%               physics parameters (e.g. GAMMA, PSPECIES_N*) to gradually
%               move from the fluid case to the kinetic case. 
dataRwmKinetic = [
%CASE# CFBAL  BETAN  q0     qa     NW  TALPHA1(Re:Im) gamma*tauA(Re:Im)  
 1     3.3218 4.1628 1.2997 7.6750 21  6.5e-4 2.0e-5  6.81963E-04  1.31584E-04
];

% RmpLin: compute plasma response to RMP fields
% Consider the base equilibrium with CFBAL=1.0, and compute linear
% steady state response of the (stable) plasma.
% Consider two sets of coils (upper and lower w.r.t. midplane) located
% just inside the wall, with artificial coil geometry.  
% Consider static RMP fields.
% input file: PROFDEN_EXAMPLE.IN, PROFROT_EXAMPLE.IN, 
%             OUTRMAR, OUTVMAR
% Here we also demonstrate how to run MARS-* in a coordinate system 
% with geometric angle, but with direct output of the amplitude of
% resonant harmonics at rational surfaces in the SFL coordinates.
% Step 1: run CHEASE with SFL coordinates (NEGP=0,NER=2)
% Step 2: run <Rmars_EXAMPLE_RmpLin_temp> with NCONVB1=1,NCONVCS= 0.
%         save PANGLE.OUT to Data/PANGLE_EXAMPLE_PEST  
% Step 3: run CHEASE with equal-arc coordinates (NEGP=-1,NER=1)
% Step 4: run <Rmars_EXAMPLE_RmpLin_temp> with NCONVB1=1,NCONVCS= 1.
%         save PANGLE.OUT to Data/PANGLE_EXAMPLE_GEOM  
% Step 5: run the following cases (vacuum or plasma response) with 
%         NCONVB1=2,NCONVCS= 1, and extract the amplitude of resonant
%         harmonics from log_mars (search for "OUTPUT TPSIM").
% CASE#1: run <Rmars_EXAMPLE_RmpLiu_vac> to compute vacuum field
%         output: log_mars_vac + dataRmpLin_vac (as example)
% CASE#2: run <Rmars_EXAMPLE_RmpLiu_pls> to compute plasma response
%         consider single fluid resistive response (as example)
%         output: log_mars_pls + dataRmpLin_pls + dataRmpLin_torq (as example)
%                 The net torques are found from log_mars_pls by 
%                 searching for "TORQUE". The two values for the JXB
%                 torque are computed in two different ways.  
dataRmpLin_vac = [
%m          s_res          Re(b1res)      Im(b1res)
 0.20E+01   0.622020E+00   0.594732E-12  -0.551838E-04
 0.30E+01   0.849908E+00  -0.992676E-13  -0.985942E-04
 0.40E+01   0.933752E+00   0.136046E-11  -0.125484E-03
 0.50E+01   0.974014E+00  -0.848418E-12  -0.141731E-03
 0.60E+01   0.995813E+00   0.105872E-11  -0.152438E-03
];
dataRmpLin_pls = [
%m          s_res          Re(b1res)      Im(b1res)
 0.20E+01   0.622020E+00   0.118843E-07   0.676565E-07
 0.30E+01   0.849908E+00   0.179762E-05   0.130732E-05
 0.40E+01   0.933752E+00   0.884428E-05   0.871991E-06
 0.50E+01   0.974014E+00   0.459263E-04  -0.165476E-04
 0.60E+01   0.995813E+00   0.504397E-04   0.280554E-04
];
dataRmpLin_torq = [
%JXB(Nm)       NTV(Nm)       REYNOLDS(Nm)
-0.559199E+01 -0.884841E+00 -0.536306E-01
];

% RmpQLin: modeling RMP induced rotation braking
% Consider the same equilibrium and coil configuration as in RmpLin.
% input file: PROFDEN_EXAMPLE.IN, PROFROT_EXAMPLE.IN, 
%             PANGLE_EXAMPLE_PEST, PANGLE_EXAMPLE_GEOM,
%             OUTRMAR, OUTVMAR
% run: Rmars_EXAMPLE_RmpQLin
%      include resonant EM torque(JXB),NTV torque and Reynolds stress (REY)
%      with Dirichlet boundary condition for the momentum equation (ITSATURAT=2)
%      note for this sample run, the computation is terminated after 20 time steps.
% output: TIMEEVOL_sample
%         The basic time evolution data are stored in the MARS-Q outputfile
%         <TIMEEVOL.OUT>.


% Tm: computing growth rate of TM
% run Rchease_EXAMPLE_TM (increase radial resolution near rational surfaces) 
% As an example, consider a resistive plasma without flow, and the
% TM in the presence of an ideal wall located at rw=1.2308a (NV=NW-1)
% Compute the TM stability in a SFL coordinate system.
% As an example, consider a uniform resistivity
% input file: PROFDEN_EXAMPLE.IN, 
%             OUTRMAR, OUTVMAR
% run: Rmars_EXAMPLE_Tm
% CASE#1:  zero beta plasma: CFBAL=0.0
%          scan ETA=1/S
% CASE#2:  fix ETA=1.0e-6, scan CFBAL
dataTm_1 = [
%CFBAL GEXP q0     qa     NW ETA     TALPHA1(Re:Im) gamma*tauA(Re:Im)
 0.0   0.0  1.5728 5.8772 18 1.0e-5  5.0e-3 0.0e-0  2.05610E-03  8.53587E-14
 0.0   0.0  1.5728 5.8772 18 1.0e-6  7.0e-4 0.0e-0  6.56848E-04 -4.12334E-14
 0.0   0.0  1.5728 5.8772 18 1.0e-7  2.5e-4 0.0e-0  1.86420E-04 -1.83179E-14
 0.0   0.0  1.5728 5.8772 18 1.0e-8  0.5e-4 0.0e-0  4.94243E-05 -4.79890E-15
];
data = dataTm_1;
x    = 1./data(:,6);
y    = data(:,9);
yy   = x.^(-3/5);
yy   = yy/yy(end)*y(end);
figure(3)
loglog(x,y,'ro',x,yy,'b--')

dataTm_2 = [
%CFBAL GEXP   q0     qa     NW ETA     TALPHA1(Re:Im) gamma*tauA(Re:Im)
 0.00  0.0000 1.5728 5.8772 18 1.0e-6  7.0e-4 0.0e-0  6.56848E-04 -4.12334E-14
 0.01  0.0120 1.5718 5.8803 18 1.0e-6  6.5e-4 0.0e-0  4.00512E-04 -9.68959E-15
 0.05  0.0599 1.5678 5.8928 18 1.0e-6  6.5e-4 0.0e-0  2.99661E-04  2.83656E-15
 0.1   0.1201 1.5629 5.9088 18 1.0e-6  2.0e-4 0.0e-0  1.67563E-04  8.74128E-15
 0.11  0.1321 1.5619 5.9120 18 1.0e-6  1.6e-4 0.0e-4  1.30830E-04  3.21197E-14
 0.115 0.1381 1.5614 5.9137 18 1.0e-6  1.3e-4 0.0e-4  1.04579E-04 -1.45784E-14
 0.12  0.1441 1.5609 5.9153 18 1.0e-6  0.8e-4 0.3e-4  7.08363E-05  2.89012E-05
 0.13  0.1562 1.5600 5.9185 18 1.0e-6  0.7e-4 0.5e-4  6.02830E-05  6.48705E-05
 0.14  0.1683 1.5590 5.9218 17 1.0e-6  0.5e-4 0.8e-4  5.01016E-05  8.73176E-05
 0.15  0.1803 1.5580 5.9251 17 1.0e-6  0.5e-4 1.0e-4  4.03767E-05  1.05733E-04
 0.17  0.2045 1.5561 5.9317 17 1.0e-6  0.2e-4 1.4e-4  2.07981E-05  1.37369E-04
 0.19  0.2287 1.5541 5.9383 17 1.0e-6  1.0e-6 1.6e-4  2.90348E-06  1.65804E-04
 0.2   0.2408 1.5531 5.9417 17 1.0e-6 -1.0e-5 1.8e-4 -5.15639E-06  1.79244E-04
 0.22  0.2650 1.5512 5.9484 17 1.0e-6 -1.5e-5 2.1e-4 -1.98477E-05  2.04450E-04
 0.24  0.2893 1.5493 5.9551 17 1.0e-6 -4.0e-5 2.3e-4 -3.41172E-05  2.27305E-04
 0.26  0.3135 1.5474 5.9619 17 1.0e-6 -5.0e-5 2.5e-4 -4.82595E-05  2.48515E-04	
 0.28  0.3378 1.5454 5.9688 17 1.0e-6 -5.5e-5 2.7e-4 -6.26146E-05  2.71614E-04
 0.3   0.3621 1.5435 5.9757 17 1.0e-6 -7.5e-5 3.0e-4 -8.50781E-05  2.97392E-04
];
data = dataTm_2;
x    = data(:,9);
y    = data(:,10);
figure(4)
plot(x,y,'b-o')

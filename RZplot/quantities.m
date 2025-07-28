% basic constants (SI)
e = 1.6021917e-19;  %[C]
me = 9.1095e-31;    %[kg]
mi = 1.67261e-27;   %[kg]
epsilon0 = 8.8542e-12; %[C^2/(N m^2)]
mu0 = pi*4e-7;      %[Wb/(A m)]
c = 2.99792e+8; %[m/s]
T0 = -273.15; %[K], absolute temperature
kB = 1.38e-23; %[J/K], Boltzmann constant
lnL = 15.0;  

if 1==0
%SDIR = '/cscratch/liuy/WorkIR_D/CU/';
%SDIR = '/cscratch/liuy/WorkD3D184003/';
%SDIR = '/home/liuy/Work/EAST/94088ke/';
SDIR = '/home/liuy/Work/DIII-D/179358/2D/';
data = load([SDIR 'global.txt']);

% AUG
aug.name = 'DIII-D';
aug.A = data(4);
aug.R = data(1);
aug.B = data(2);  
aug.Te = data(6); %[eV]
aug.Ti = data(7);
aug.n  = data(5);
aug.Z  = 2;
aug.Zeff = 1.0; 

d = load([SDIR '../PROFROT_SAVE']);
rot0 = d(2,2);

else

aug.name = 'D3D';
aug.A = 2.8524;
aug.R = 1.6955;
aug.B = 2.0874;  
aug.Te = 3.4383e+3; %[eV]
aug.Ti = 7.0184e+3;
aug.n  = 7.2886e+19;
aug.Z  = 2;
aug.Zeff = 1.0; 
rot0   = 1.2897e+5;
end


% calculate all quantities
dev = aug;
name = dev.name;
A = dev.A;
R = dev.R;
B = dev.B;
Te = dev.Te;
Ti = dev.Ti;
n = dev.n;
Z = dev.Z;
Zeff = dev.Zeff;
a = R/A;
Te = Te*e/kB; %[K]
Ti = Ti*e/kB; %[K]
vthi = sqrt(2*kB*Ti/mi);
vthe = sqrt(2*kB*Te/me);

disp(['quantities for ', name])
% gyrofrequencies
omega_ce = e*B/me; 
omega_ci = e*B/mi/Z; %[Hz]
disp(['f_ci=', num2str(omega_ci/2/pi), ' f_ce=', num2str(omega_ce/2/pi)]), 

% plasma frequency
omega_pi = sqrt(e^2*n/mi/epsilon0);
omega_pe = sqrt(e^2*n/me/epsilon0);
disp(['f_pi=',num2str(omega_pi/2/pi), ' f_pe=',num2str(omega_pe/2/pi)])

% toroidal gyro-radius
r_Le = vthe/(omega_ce); %[m]
r_Li = vthi/(omega_ci); %[m]
disp(['r_Li=', num2str(r_Li), ' r_Le=', num2str(r_Le)]), 

% poloidal gyro-radius/a, assuming q=3
q=3;  
r_Le_p = r_Le*q*A/a; %[m]
r_Li_p = r_Li*q*A/a; %[m]
disp(['r_Li_p/a=', num2str(r_Li_p), ' r_Le_p/a=', num2str(r_Le_p)]), 

% Debye length
r_D = sqrt(epsilon0*kB*Te/n/e^2); %[m]
disp(['r_D=', num2str(r_D)]), 

% rho_*
rho_star = r_Li/a;
% k_perp
k_perp = [1 15]*100; %[1/m]
disp(['k_perp=', num2str(k_perp), ' rho_*=', num2str(rho_star), ' k_perp*rho_*=' num2str(k_perp*rho_star)]), 

% Alfven frequency
v_A = B/sqrt(mu0*mi*Z*n); %[m/s]
omega_A = v_A/R;  
tau_A = R/v_A; %[s]
disp(['f_A=',num2str(omega_A/2/pi), ' tau_A=',num2str(tau_A)])

% vthi/v_A
disp(['vthi/v_A=',num2str(vthi/v_A)])

OMEGACI0 = omega_ci/omega_A

% ion sound frequency
Gamma = 5/3;
rho = mi*Z*n;
P = n*kB*(Ti+Te);
v_s = sqrt(Gamma*P/rho);
omega_s = v_s/R;
tau_s = R/v_s;
disp(['f_s=',num2str(omega_s/2/pi), ' tau_s=',num2str(tau_s)])

% diamagnetic drift frequency, estimation!
omega_stari = 2*n*kB*Ti/a/(e*n*B)/R;
omega_stare = 2*n*kB*Te/a/(e*n*B)/R;
disp(['f_*i=',num2str(omega_stari/2/pi),' f_*e=',num2str(omega_stare/2/pi)])

% magnetic drift frequency
q = 2;
omega_Bi = q*vthi^2/omega_ci/R/a;
omega_Be = q*vthe^2/omega_ce/R/a;
disp(['f_Bi=',num2str(omega_Bi/2/pi), ' f_Be=',num2str(omega_Be/2/pi)])

% transit frequency (circulating ions and electrons)
omega_ti = vthi/q/R;
omega_te = vthe/q/R;
disp(['f_ti=',num2str(omega_ti/2/pi), ' f_te=',num2str(omega_te/2/pi)])

% bounce frequency (trapped ions and electrons)
omega_bi = omega_ti/sqrt(2*A);
omega_be = omega_te/sqrt(2*A);
disp(['f_bi=',num2str(omega_bi/2/pi), ' f_be=',num2str(omega_be/2/pi)])

%(thermal) particle collision 
nu_e = 2.91e-12*n*lnL/(Te*kB/e)^1.5; %collision rate in [1/s]
nu_i = 4.80e-14*Zeff^4/sqrt(Z)*n*lnL/(Ti*kB/e)^1.5;  %[1/s]
ct_e = 1/nu_e;  %collision time gap [s]
ct_i = 1/nu_i;  %collision time gap [s]
mfp_e = vthe*ct_e;  %mean-free-path [m]
mfp_i = vthi*ct_i;  %mean-free-path [m]
disp(['collision rate [1/s]: nu_e=' num2str(nu_e) ', nu_i=',num2str(nu_i)])
disp(['collision time gap [s]: ct_e=' num2str(ct_e) ', ct_i=',num2str(ct_i)])
disp(['mean-free-path [m]: mfp_e=' num2str(mfp_e) ', mfp_i=',num2str(mfp_i)])

%Izzo term for collisionality
%grel = (1-(vthe/c)^2)^(-0.5);
%f_Izzo = e^3*lnL*n/(grel*vthe^2*4*pi*epsilon0^2*me)

% plasma resistivity according to Spitzer's formula
%C0   = 4*sqrt(2*pi)/3; %=3.3422
C0   = pi;
etaP = C0*Zeff*e^2*sqrt(me)*lnL/(4*pi*epsilon0)^2/(kB*Te)^1.5; %eta_perp
eta0 = mu0*R^2*omega_A/A^2;
S    = eta0/etaP;
ETA  = 1/S;
Ic   = 1e+3/(R*B/4e-7/pi);
ROTE = rot0/omega_A;
disp(['eta=',num2str(etaP),'[Ohm m];  S=',num2str(S)])
format short e
RES = [tau_A S ETA ROTE Ic]

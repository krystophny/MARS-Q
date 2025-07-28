% comput Va, Vti, rot, and nui
function [Vti,nui,rot,peqi] = MacGetNTVV(Rmid)

global Mac
global SDIR

%constants
e = 1.6021917e-19;  %[C]
mi = 1.67261e-27;   %[kg]
mu0 = pi*4e-7;      %[Wb/(A m)]
eps0 = 8.8542e-12;  %[F/m] 
kB = 1.38e-23;      %[J/K], Boltzmann's constant

if Mac.dataEXP>0
  eval(['load ' SDIR 'dataDEN'])
  eval(['load ' SDIR 'dataTEM'])
  eval(['load ' SDIR 'dataROT'])
end

%JET plasma parameters
B0 = Mac.B0EXP;
R0 = Mac.R0EXP;
ne = 2.9e+19;  if Mac.dataEXP>0, ne=dataDEN(1,2); end
Zi = 1;
mui= 2;
Zeff = 2.865;
ni = (12-Zeff)/(5*Zeff)*ne;
mi  = mi*mui;

%Alfven velocity
rho = mi*ni;
v_A = B0/sqrt(mu0*rho);  %Alfven velocity in [m/s]
f_A = v_A/R0/2/pi;

%pressure and density profiles
if Mac.dataEXP==0 
  peq = Mac.peq;
  den = Mac.den; 
  rot = Mac.rot;
  peqi= peq/2;
else
  Ti = (dataTEM(:,Mac.dataEXP+1) + dataTEM(:,Mac.dataEXP+2))/2;
  den = dataDEN(:,2)/dataDEN(1,2);
  peq = 2*ni*den.*Ti*e*mu0/B0^2;
  peqi= peq/2;
  rot = (dataROT(:,Mac.dataEXP+1) + dataROT(:,Mac.dataEXP+2))/2/(v_A/R0);
end

%test with Tim [nu_i(0)=1.0584e+3Hz]
%ni = 3.0e+19*ones(size(ni));
%Ti = 5.0e+3;
%peq = 2*ni*Ti*e*mu0/B0^2;

%ion thermal velocity Vti
if Mac.dataEXP == 0
  Vti = sqrt(peq./den);
else
  Vti = sqrt(2*Ti*e/mi)/v_A;
end

if Mac.plot_Vti
  figure(Mac.plot_Vti)
  plot(Rmid,Vti,'b-','LineWidth',2), hold on,
  xlabel('R [m]','FontSize',16),
  ylabel('ion thermal velocity/v_A','FontSize',14),
end

%lnLamb
%lnLamb = 20.39 - log(2*kB^1.5*B0/sqrt(mu0)/v_A^3/mi^1.5/Zi^1.5*den.^1.5./peq); %CHEASE formula
lnLamb = 23-log(4.0e-3*B0*Zi^3*e^1.5*den.^2/mu0^0.5/mi^2/v_A^4./peq.^1.5); %NRL formulary
lnLamb(end) = lnLamb(end-1);

%ion-collisionality nui
%nui = 4.0*pi*B0^2*kB^1.5*e/mi^3/v_A^5/mu0/Zi^1.5*den.^2.5.*lnLamb./peq.^1.5; %Shaing
%nui = 4.8e-14*2^1.5*Zi^4/mui^0.5/mu0*e^1.5*B0^2/mi^2.5/v_A^5.*den.^2.5./peq.^1.5.*lnLamb; %NRL
nui = 1.0/eps0^2/mi^3*Zi^4*e^4*B0^2/(4*pi*mu0*v_A^5)*lnLamb.*den.^2.5./peq.^1.5; %Sabbagh
nui(1)
nui = nui/f_A;
nui(end) = nui(end-1);

%re-normalization
Vti  = Vti*v_A;
nui  = nui*f_A;
rot  = rot*f_A;
peqi = peqi*Mac.B0EXP^2/mu0;

%ion collisionality
nu_star = sqrt(2)*sqrt(pi)/3*Mac.q.*Rmid.*nui./Vti;

if Mac.plot_nui > 0
  figure(Mac.plot_nui)
  subplot(2,1,1), plot(Mac.s(1:Mac.Ns1),lnLamb,'LineWidth',2), hold on,
                  ylabel('ln\Lambda','FontSize',16),
  subplot(2,1,2), plot(Mac.s(1:Mac.Ns1),nui,'LineWidth',2), hold on,
                  xlabel('s','FontSize',16),
                  ylabel('\nu_i','FontSize',16),
end


if Mac.plot_NTVCnu > 0
  figure(Mac.plot_NTVCnu)
  plot(Rmid,nu_star,'b-','LineWidth',2), hold on,
  xlabel('R [m]','FontSize',16),
  ylabel('C_\nu','FontSize',16),
end

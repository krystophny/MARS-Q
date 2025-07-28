% Read MARS-Q output data from the NTV module
%   PROFNTV.OUT   :  relevant quantities for NTV computation  
%   TORQUENTV     :  NTV torque density
% Plot the above information
% Note that this is a stand-alone script w/o calling other Mac-procedures

% read in data
if kplot_tntv(1)>0, torqn= load([SDIR 'TORQUENTV.OUT']);  end
if kplot_tntv(2)>0, torq = load([SDIR 'TORQUENTV2.OUT']); end 
if kplot_tntv(3)>0, prof = load([SDIR 'PROFNTV.OUT']);    end
if kplot_tntv(4)>0, rot  = load([SDIR 'PROFOFFSET.OUT']); end

fT   = 1.0;
fV   = 1.0;
if kscale==1
   fT = B0^2/(4e-7*pi);  %10 kAt ==> factor of 100 here
   fV = 4*pi^2*R0^3/ASPECT^2;
end


if kplot_tntv(3) > 0
% profiles
rho   = prof(:,1);  %radial coordinate
q     = prof(:,2);  %safety factor
eps   = prof(:,3);  %equivalent r/R0 in Shaing's theory
rhoi  = prof(:,4);  %thermal ion density
nu    = prof(:,5);  %collisionality, thermal ion or electron depending on ZCHARGE
wti   = prof(:,6);  %thermal ion transit frequency
wB0   = prof(:,7);  %gradB precession drift frequency, depending on ZCHARGE
qwe   = prof(:,8);  %ExB drift frequency, =q*omega_E
qws   = prof(:,9);  %diamagnetic drift frequency, =q*omega_*p, depending on ZCHARGE 
qwsT  = prof(:,10); %grad T drift frequency, =q*omega_*T, depending on ZCHARGE
dB    = prof(:,11); %surface averaged |dB|/B0

psi   = rho.^2;

% some fix-ups
dB(1) = dB(3);
dB(2) = dB(3);

% calculate "boundary"-frequencies between different NTV regimes
% first for non-resonant NTV torque
%           nu < nu_n1 : nu-regime
%   nu_n1 < nu < nu_n2 : sqrt(nu)-regime
%   nu_n2 < nu < nu_n3 : 1/nu-regime
%   nu_n3 < nu         : violation of assumptions in Shaing's NTV theory         
nu_n1 = abs(qwe./q).*(dB./eps).^2;
nu_n2 = abs(qwe./q);
nu_n3 = sqrt(eps).*wti; 

% next for resonant NTV torque
%           nu < nu_r1 : superbanana-regime
%   nu_r1 < nu < nu_r2 : superbanana-plateau-regime
%   nu_r2 < nu < nu_r3 : 1/nu-regime
%   nu_r3 < nu         : violation of assumptions in Shaing's NTV theory         
nu_r1 = wB0.*(dB./eps).^1.5;
nu_r2 = wB0./eps;
nu_r3 = sqrt(eps).*wti;

% plot the equivalent minor radius eps used in Shaing's NTV calculations
hf = figure(10*kplot_tntv(3)+1);
plot(rho,eps,'b-','LineWidth',2), hold on,
xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
ylabel('\epsilon_{NTV}','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)

% plot ExB frequency
hf = figure(10*kplot_tntv(3)+2);
plot(psi,qwe./q,'b-','LineWidth',2), hold on,
xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
ylabel('ExB frequency [rad/s]','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
%legend('\omega_E^{EXPT.}','\omega_E^{MARS-Q}')

res = [rho qwe./q*TAUA];
save PROFWE res -ascii -double

% plot boundaries for various NTV regimes
% together with the collision frequency nu
hf = figure(10*kplot_tntv(3)+3);
semilogy(psi,nu_n1,'b:','LineWidth',1), hold on,
semilogy(psi,nu_n2,'b--','LineWidth',2), hold on,
semilogy(psi,nu_n3,'b-','LineWidth',3), hold on,
semilogy(psi,nu_r1,'r:','LineWidth',1), hold on,
semilogy(psi,nu_r2,'r--','LineWidth',2), hold on,
semilogy(psi,nu_r3,'r-','LineWidth',3), hold on,
semilogy(psi,nu./eps,'k-','LineWidth',4), hold on,
xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
ylabel('frequencies [rad/s]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
legend('n1','n2','n3','r1','r2','r3','\nu_{eff}')

% save effective collisionality for use in MARS-K
%nu_eff  = nu./eps*TAUA;
nu_eff  = nu*TAUA;  %this choice benchmarked by ZR Wang
s_new   = linspace(0,1,101)';
nu_new  = spline(sqrt(psi),nu_eff,s_new);
res     = [s_new nu_new];
save PROFNUI res -ascii

hf = figure(kplot_tntv(3)+4);
plot(psi,dB,'b-','LineWidth',3), hold on,
xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
ylabel('<|\delta B|>/B_0','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)

end

if kplot_tntv(2) > 0
% get NTV torque densities
s    = torq(:,1);  %radial mesh
Ttot = torq(:,2);  %total NTV torque density = Tnon + Tres
Tnon = torq(:,3);  %non-resonant portion of NTV torque density
Tres = torq(:,4);  %resonant portion of NTV torque density
psi  = s.^2;

% plot NTV torque densities
hf = figure(kplot_tntv(2));
res_ntv = [psi (Ttot)*fT];
if kaxisy==1
   plot(psi,(Tnon)*fT,'b-','LineWidth',2), hold on,
   plot(psi,(Tres)*fT,'r-','LineWidth',2), hold on,
   plot(psi,(Ttot)*fT,'k-','LineWidth',3), hold on,
else
   semilogy(psi,abs(Tnon)*fT,'b--','LineWidth',2), hold on,
   semilogy(psi,abs(Tres)*fT,'k--','LineWidth',2), hold on,
   semilogy(psi,abs(Ttot)*fT,'r-','LineWidth',3), hold on,
end
xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
if kscale==1
   ylabel('NTV torque density [N/m^2]','FontSize',16,'FontWeight','Bold')
else
   ylabel('normalised NTV torque density','FontSize',16,'FontWeight','Bold')
end
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
%legend('T_{JXB}','T_{NTV}^{non}','T_{NTV}^{res}','T_{NTV}^{non}+T_{NTV}^{res}')
legend('T_{NTV}^{non}','T_{NTV}^{res}','T_{NTV}^{non}+T_{NTV}^{res}')

end

if kplot_tntv(1) > 0

s    = torqn(:,1);
Tntv = torqn(:,2);
Dntv = torqn(:,3);
TntvI= torqn(:,6);
TntvE= torqn(:,7);

if ksmooth==1
   s0 = sqrt(psipp(ksmooth));
   ff = 1-0.95*exp(-(s-s0).^2/0.001^2);
   Tntv = Tntv.*ff;
   TntvI = TntvI.*ff;
   TntvE = TntvE.*ff;
elseif ksmooth==2
   c = 0.2;
   for k=1:20
       Tntv(2:end-1) = Tntv(2:end-1)*(1-2*c) + Tntv(1:end-2)*c + Tntv(3:end)*c; 
   end
end

if kaxisx==3
   d = load([SDIR 'PROFEQ.OUT']);
   q = d(:,2);
   x = (q(1:end-1)+q(2:end))/2;
elseif kaxisx==2
   x = s.^2;
elseif kaxisx==1
   x = s;
end

hf = figure(kplot_tntv(1));
if kaxisy==1
   %plot(x,Tntv*fT,'r-',x,TntvI*fT,'b-',x,TntvE*fT,'k-','LineWidth',3), hold on,
   plot(x,Tntv*fT,'r-','LineWidth',3), hold on,
else
   %semilogy(x,abs(Tntv)*fT,'r-',x,abs(TntvI)*fT,'b-',x,abs(TntvE)*fT,'k-','LineWidth',3), hold on,
   %semilogy(x,-TntvI*fT,'r-',x,TntvE*fT,'b-','LineWidth',3), hold on,
   semilogy(x,abs(Tntv)*fT,'b:','LineWidth',3), hold on,
end
if 1==0
   a=axis; 
   for k=1:length(psipp)
       if kaxisy==1, plot([psipp(k) psipp(k)],[a(3) a(4)],'k--'), end
       if kaxisy==2, semilogy([psipp(k) psipp(k)],[a(3) a(4)],'k--'), end
   end
   format short e
   d  = load([SDIR 'PROFEQ.OUT']);
   q  = d(:,2);
   qx = d(:,1).^2;
   psipq = spline(qx,q,psipp)
end

if kaxisx==3
   xlabel('q','FontSize',18,'FontWeight','Bold')
elseif kaxisx==2
   xlabel('\psi_n','FontSize',18,'FontWeight','Bold')
elseif kaxisx==1
   xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
end
if kscale==1
   if kaxisy==1
      ylabel('NTV torque [N/m^2]','FontSize',16,'FontWeight','Bold')
   else
      %ylabel('-T_{NTV} [N/m^2]','FontSize',16,'FontWeight','Bold')
      ylabel('NTV torque [N/m^2]','FontSize',16,'FontWeight','Bold')
   end
else
   ylabel('normalised NTV torque density','FontSize',16,'FontWeight','Bold')
end
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
%legend('ion+electron','ion','electron')
%legend('-T^{ion}_{NTV}','T^{ele}_{NTV}')

% calculate net NTV torque in [Nm]
NTntv  = sum( (Tntv(1:end-1)+Tntv(2:end)).*diff(s) )/2*fT*fV
NTntvI = sum( (TntvI(1:end-1)+TntvI(2:end)).*diff(s) )/2*fT*fV
NTntvE = sum( (TntvE(1:end-1)+TntvE(2:end)).*diff(s) )/2*fT*fV

%plot density pumpout fraction associated with NTV
hf = figure(10*kplot_tntv(1)+1);
if kaxisy==1
   plot(x,Dntv,'b-','LineWidth',3), hold on,
else
   semilogy(x,abs(Dntv),'b-','LineWidth',3), hold on,
end
if kaxisx==3
   xlabel('q','FontSize',18,'FontWeight','Bold')
elseif kaxisx==2
   xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
elseif kaxisx==1
   xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
end
ylabel('density pumpout due to NTV','FontSize',14,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)

end

if kplot_tntv(4) > 0
% plot offset rotation profile
s  = rot(:,1);
w1 = rot(:,2);  %offset flow
w0 = rot(:,3);  %initial flow

hf = figure(kplot_tntv(4));
plot(s.^2,w0/TAUA,'r--',s.^2,w1/TAUA,'b-','LineWidth',3), hold on,
xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
ylabel('\Omega/\omega_A','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
legend('initial flow','offset flow')

% generate PROFROT file using the offset flow
ROTE  = max(abs(w1))

fid = fopen([SDIR 'PROFROT_NEW'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; w1']);
fclose(fid);

% plot omegaE from D3D data
if kdev==3 & 1==0
data = load([SDIR 'PROFWE']);
wE   = spline(data(:,1),data(:,2),s);
figure(2)
plot(s.^2,wE,'r-','LineWidth',3)
end

end

% plot torque density from MARS-K computed energy
if kdev==12 & kaxisx==2 & kaxisy==2 & 1==1
data = load([SDIR 'ENERGY_DENSITY']);
xx   = data(:,1).^2;
yy   =-data(:,3)*2;
figure(kplot_tntv(1))
semilogy(xx,abs(yy),'b--','LineWidth',3)
end
 

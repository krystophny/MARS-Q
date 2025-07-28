% Read MARS-Q output data for plotting the REY torque density
%   TORQUEREY     :  Reynolds stress torque density
% Note that this is a stand-alone script w/o calling other Mac-procedures

% read in data
torq = load([SDIR 'TORQUEREY.OUT']);

% get REY torque densities
s    = torq(:,1);  %radial mesh
Trey = torq(:,2);  %REY torque density 
Drey = torq(:,3);  %density pumpout fraction due to MHD

fT   = 1.0;
if kscale==1
   fT = B0^2/(4e-7*pi);
elseif kscale==2
   fT = 1/max(abs(Trey));
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

% plot REY torque densities
hf = figure(kplot_trey);

res_rey = [x (Trey)*fT];

if kaxisy==1
   plot(x,(Trey)*fT,'k-','LineWidth',3), hold on,
else
   semilogy(x,abs(Trey)*fT,'k-','LineWidth',3), hold on,
end
if kaxisx==3
   xlabel('q','FontSize',18,'FontWeight','Bold')
elseif kaxisx==2
   xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
elseif kaxisx==1
   xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
end
if kscale==1
   ylabel('Reynolds stress torque density [N/m^2]','FontSize',16,'FontWeight','Bold')
else
   ylabel('normalised Reynolds stress torque density','FontSize',16,'FontWeight','Bold')
end
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

% plot density pumpout fraction due to MHD
hf = figure(10*kplot_trey+1);

if kaxisy==1
   plot(x,(Trey),'k-','LineWidth',3), hold on,
else
   semilogy(x,abs(Trey),'k-','LineWidth',3), hold on,
end
if kaxisx==3
   xlabel('q','FontSize',18,'FontWeight','Bold')
elseif kaxisx==2
   xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
elseif kaxisx==1
   xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
end
ylabel('density pumpout due to MHD','FontSize',14,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')



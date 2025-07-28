% Read MARS-Q output data for plotting the JXB torque density
%   TORQUEJXB     :  JXB torque density
% Note that this is a stand-alone script w/o calling other Mac-procedures

% read in data
torq = load([SDIR 'TORQUEJXB.OUT']);

% get NTV torque densities
s    = torq(:,1);  %radial mesh
Tjxb = torq(:,2);  %JXB torque density 

fT   = 1.0;
if kscale==1
   fT = B0^2/(4e-7*pi);
elseif kscale==2
   fT = 1/max(abs(Tjxb));
end

% plot JXB torque densities
% note that the torque density is normalised by fac=B0^2/mu0
hf = figure(kplot_tjxb);

if kaxisx==3
   d = load([SDIR 'PROFEQ.OUT']);
   q = d(:,2);
   x = (q(1:end-1)+q(2:end))/2;
elseif kaxisx==2
   x = s.^2;
elseif kaxisx==1
   x = s;
end

res_jxb = [x (Tjxb)*fT];

if kaxisy==1
   plot(x,(Tjxb)*fT,'r-','LineWidth',3), hold on,
else
   semilogy(x,abs(Tjxb)*fT,'r-','LineWidth',3), hold on,
end
if kaxisx==3
   xlabel('q','FontSize',18,'FontWeight','Bold')
elseif kaxisx==2
   xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
elseif kaxisx==1
   xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
end
if kscale==1
   ylabel('JXB torque density [N/m^2]','FontSize',16,'FontWeight','Bold')
else
   ylabel('normalised JXB torque density','FontSize',16,'FontWeight','Bold')
end
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')


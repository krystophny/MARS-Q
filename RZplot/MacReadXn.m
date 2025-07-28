% Read normal displacement of plasma surface from MARS-Q output
% PROFDISP.OUT  
% Note that this is a stand-alone script w/o calling other Mac-procedures

% read in data
dispn = load([SDIR 'PROFDISP.OUT']);

% get NTV torque densities
th  = dispn(:,1);  %poloidal mesh
xnr = dispn(:,2);  %real part of Xn
xni = dispn(:,3);  %imaginary part of Xn

% convert from [0,2*pi] mesh to [-pi,pi] mesh
II     = find(th>pi);
th(II) = th(II) - 2*pi;
[th,II] = sort(th);
xnr     = xnr(II);
xni     = xni(II);

fT   = 1.0;
if kscale==1
   fT = R0*1000;
end

% plot normal displacement
hf = figure(kplot_disp);
plot(th*180/pi,abs(xnr+xni*i)*fT,'b-','LineWidth',3), hold on,
xlabel('equal-arc poloidal angle [degree]','FontSize',16,'FontWeight','Bold')
if kscale==1
   ylabel('|normal displacement| [mm]','FontSize',16,'FontWeight','Bold')
else
   ylabel('|normal displacement|','FontSize',16,'FontWeight','Bold')
end
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)


% plot quantities from MARS-Q output
% also plot density pumpout fractions
format short e

kscale        = 0;   %=0: plot dimensionless torque; =1: plot in physical units
kaxisx        = 1;   %=1: x=s; =2: x=s^2=psip; =3: x=q & need PROFEQ.OUT
kaxisy        = 2;   %=1: linear-scale for y-axis; =2: log-scale for y-axis 
ksmooth       = 0;   %=1: smooth the peak number ksmooth from psipp-array
                     %=2: generic smoothing along radial coordinate
kplot_cont    = 0;   %need PROFEQ.OUT     if >0
kplot_tjxb    = 1;   %need TORQUEJXB.OUT  if >0  plot JXB torque density
kplot_tntv(1) = 2;   %need TORQUENTV.OUT  if >0, plot NTV torque density and density pumpout fraction due to NTV
kplot_tntv(2) = 0;   %need TORQUENTV2.OUT if >0, plot resonant vs. non-resonant NTV torque density
kplot_tntv(3) = 0;   %need PROFNTV.OUT    if >0, plot detailed profiles associated with NTV torque
kplot_tntv(4) = 0;   %need PROFOFFSET.OUT if >0, plot offset rotation profile due to NTV
kplot_trey    = 3;   %need TORQUEREY.OUT  if >0, plot torque density and density pumpout fraction due to perturbed velocity
kplot_disp    = 0;   %need PROFDISP.OUT   if >0, plot normal displacement of the plasma surface
 
kdev = 67;
MacReadParam

if kplot_cont > 0,      MacReadContinua, end
if kplot_tjxb > 0,      MacReadJXB,      end
if sum(kplot_tntv) > 0, MacReadNTV,      end
if kplot_trey > 0,      MacReadREY,      end
if kplot_disp > 0,      MacReadXn,       end


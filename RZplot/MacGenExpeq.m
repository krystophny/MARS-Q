% generate EXPEQ file with analytic profiles and plasma boundary shape
% 2009-09-25, YQ Liu

% *****************
% input parameters
% *****************
% input parameters for specifying plasma and wall boundarys
R0 = 3.0;     %aspect ratio
EL = 1.6;     %elongation
DL = 0.3;     %triangularity
XA = 0.00;    %X-point shaping parameter
XB = 100.0;   %X-point shaping parameter
XC = 2.0;     %X-point shaping parameter
RW = 1.3;     %conformal wall minor radius
R0EXP = 3.0;  %plasma major radius

% input parameters for plasma density profile
Nedge = 0.01; %edge density
N0    = 0.5;  %pedestal top density
s0    = 0.95; %pedestal top minor radius

% input parameters for (electron&ion) temperature profile 
Tedge = 0.01; %edge temperature
s1    = 1.2*s0;

% input parameters for plasma pressure
P0 = 0.5;    %plasma central pressure, normalised by B0EXP^2/mu0

% input parameters for plasma toroidal current density 
J0 = 2.0;     %plasma central current density, normalised by B0EXP/R0EXP/mu0
JV = 1.0;     %current shape parameter

% input parameters for toroidal rotation profile 
Wedge = 0.05; %edge rotation

% **************
% creating data
% **************
% plasma boundary and wall
t  = linspace(-pi,pi,301)';
g  = 1.0 + XA*exp(-XB*((t+pi/2)/(pi/2)).^XC);
R  = 1.0 + cos(t+DL*sin(t))/R0;
Z  = EL*sin(t)/R0.*g;
R2 = 1.0 + (R-1.0)*RW;
Z2 = Z*RW;

hf = figure;
plot(R*R0EXP,Z*R0EXP,'r-',R2*R0EXP,Z2*R0EXP,'b-','LineWidth',3), hold on,
xlabel('R [m]','FontSize',16,'FontWeight','Bold')
ylabel('Z [m]','FontSize',16,'FontWeight','Bold')
ha = get(hf,'CurrentAxes');  set(ha,'FontSize',16,'FontWeight','Bold')
axis equal

% plasma minor radius
s = linspace(0,1,201)';

% density profile with pedestal
f     = zeros(size(s));
II    = find(s>s0);
f(II) = (s(II)-s0).^2./(1-s0)^2;
N     = ( (1-(1-N0)*(s/s0).^2).*(1-f) + Nedge )/(1+Nedge);
 
hf = figure;
plot(s,N,'b-','LineWidth',3), hold on,
plot([s0 s0],[0 1],'k--'), hold on,
ylabel('normalised density','FontSize',16,'FontWeight','Bold')
xlabel('s','FontSize',16,'FontWeight','Bold')
ha = get(hf,'CurrentAxes');  set(ha,'FontSize',16,'FontWeight','Bold')

% temperature profile
T     = ( (1-s1*s.^2+s.^3/3).*(1-f) + Tedge )/(1+Tedge);

hf = figure;
plot(s,T,'b-','LineWidth',3), hold on,
plot([s0 s0],[0 1],'k--'), hold on,
ylabel('normalised temperature','FontSize',16,'FontWeight','Bold')
xlabel('s','FontSize',16,'FontWeight','Bold')
ha = get(hf,'CurrentAxes');  set(ha,'FontSize',16,'FontWeight','Bold')

% pressure profile and P'
P   = P0*N.*T;
sx  = (s(1:end-1)+s(2:end))/2;
Px  = diff(P)./diff(s)./sx/2;
PP  = spline(sx,Px,s);

hf = figure;
plot(s,P,'b-','LineWidth',3), hold on,
plot([s0 s0],[0 1],'k--'), hold on,
ylabel('normalised pressure','FontSize',16,'FontWeight','Bold')
xlabel('s','FontSize',16,'FontWeight','Bold')
ha = get(hf,'CurrentAxes');  set(ha,'FontSize',16,'FontWeight','Bold')

hf = figure;
plot(s,PP/P0,'b-','LineWidth',3), hold on,
ylabel('pressure gradient','FontSize',16,'FontWeight','Bold')
xlabel('s','FontSize',16,'FontWeight','Bold')
ha = get(hf,'CurrentAxes');  set(ha,'FontSize',16,'FontWeight','Bold')
a = axis;
plot([s0 s0],[a(3) a(4)],'k--'), hold on,

% toroidal current profile
J  = J0*(1-s.^(2*JV));

hf = figure;
plot(s,J,'b-','LineWidth',3), hold on,
ylabel('toroidal current density','FontSize',16,'FontWeight','Bold')
xlabel('s','FontSize',16,'FontWeight','Bold')
ha = get(hf,'CurrentAxes');  set(ha,'FontSize',16,'FontWeight','Bold')

% toroidal rotation profile
w = (1 - 2*s.^2 + s.^3 + Wedge)/(1+Wedge);

hf = figure;
plot(s,w,'b-','LineWidth',3), hold on,
ylabel('normalised toroidal rotation','FontSize',16,'FontWeight','Bold')
xlabel('s','FontSize',16,'FontWeight','Bold')
ha = get(hf,'CurrentAxes');  set(ha,'FontSize',16,'FontWeight','Bold')

% **************
% output data
% **************
% output density profile
data = [length(N) 1];
fid  = fopen('f1','w'); 
fprintf(fid,'%4i %4i\n',data);
fclose(fid);
data = [s N];
save f2 data -ascii -double
!cat f1 f2 > PROFDEN

% output rotation profile
data = [length(w) 1];
fid  = fopen('f1','w'); 
fprintf(fid,'%4i %4i\n',data);
fclose(fid);
data = [s w];
save f2 data -ascii -double
!cat f1 f2 > PROFROT

% output EXPEQ file
data = [1/R0 0 0]';
save f1 data -ascii -double
data = [length(R) 2 2];
fid  = fopen('f2','w'); 
fprintf(fid,'%4i %4i %4i\n',data);
fclose(fid);
data = [R Z; R2 Z2];
save f3 data -ascii -double
data = [length(s) 2]';
fid  = fopen('f4','w'); 
fprintf(fid,'%4i\n %4i\n',data);
fclose(fid);
data = [s; PP; J];
save f5 data -ascii -double
!cat f1 f2 f3 f4 f5 > EXPEQ

!rm -f f1 f2 f3 f4 f5



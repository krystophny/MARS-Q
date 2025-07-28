function B=MacACADbound(R0EXP)

global Acad

data = load('BOUNDP');

data = data*R0EXP;

N = size(data,1);

R = data(1:N,1);
Z = data(1:N,2);

Z0EXP = (min(Z)+max(Z))/2;

%plasma shape
%figure(1)
%plot(R,Z,'r-','LineWidth',2), hold on,
%axis equal
%xlabel('R [m]','FontSize',18,'FontWeight','Bold')
%ylabel('Z [m]','FontSize',18,'FontWeight','Bold')
%ha=get(1,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

%smooth plasma boundary slightly
t  = atan2(Z-Z0EXP,R-R0EXP);
[t,II] = sort(t);
R  = R(II);
Z  = Z(II);
[t,II] = unique(t);
R  = R(II);
Z  = Z(II);


if 1==0
   t1 = linspace(-pi,pi,101);
   R1 = pchip(t,R,t1);
   Z1 = pchip(t,Z,t1); 
else
   t1 = t;
   R1 = R;
   Z1 = Z;
end

if 1==0
   [Y,I] = min(Z1);
   R1(I) = 1.33;
   Z1(I) =-1.14;
end

if 1==0
   [Y,I] = max(Z1);
   R1(I) = 1.36;  Z1(I) = 0.98;
end

t2 = linspace(-pi,pi,255);
R2 = pchip(t1,R1,t2);
Z2 = pchip(t1,Z1,t2);

%smoothing plasma boundary shape
c    = 0.25;
for k=1:5
    R2(2:end-1) = (1-2*c)*R2(2:end-1) + c*(R2(1:end-2)+R2(3:end));
    Z2(2:end-1) = (1-2*c)*Z2(2:end-1) + c*(Z2(1:end-2)+Z2(3:end));
end
%plot(R1,Z1,'bx','LineWidth',2), hold on,
%plot(R2,Z2,'b-','LineWidth',2), hold on,

RP = R2'; 
ZP = Z2';

%shifting whole plasma by DRPW along R-axis
if isfield(Acad,'DRPW')
   RP = RP+Acad.DRPW*1e-2;
   R0EXP = R0EXP+Acad.DRPW*1e-2;
   Acad.R0EXP = R0EXP;
end   

%DIII-D wall surface
if length(Acad)>0
if strcmp(Acad.DEVICE,'D3D')
data = load('MacCOMMON/BOUNDW_D3D')*1.6995;
end
if strcmp(Acad.DEVICE,'MAST')
data = load('MacCOMMON/BOUNDC_MAST')*0.83222;
end
if strcmp(Acad.DEVICE,'AUG')
data = load('MacCOMMON/BOUNDC_AUG')*1.0;
end
if strcmp(Acad.DEVICE,'JET')
data = load('MacCOMMON/BOUNDW_JET')*3.0081;
end
end

RW   = data(:,1);
ZW   = data(:,2);

%plot(RW,ZW,'b'), hold on,


B = [RP ZP; RW ZW]/R0EXP;
%save BOUNDPW RZ -ascii -double


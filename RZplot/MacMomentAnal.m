% plot stead state soluton of momentum balance solution for analytic model
% from [Liu PP 20, 042503(2013)]
a0=2.0;
d0=3;
T0=3.2;
alf=1.3;
bet=2.3;
gam=1.7;

y0=T0/d0*(alf-bet)/bet/(bet+gam)*(alf/(alf+gam)-alf/(alf-bet)+gam/(alf+gam)*(exp(gam)-exp(-bet))/(exp(-alf)-exp(-bet)));

s=linspace(0,1,101);
y=y0*(alf*exp(-bet*s)-bet*exp(-alf*s))/(alf-bet)-T0/d0*alf/(bet+gam)*((exp(gam*s)-exp(-alf*s))/(alf+gam)-(exp(-bet*s)-exp(-alf*s))/(alf-bet));

plot(s,y)

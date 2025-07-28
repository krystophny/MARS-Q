function MacGetBsn(B1,B2)

global Mac
global chis Bss1 Bs10 Bs20 Bs30 chis0

[smin,II] = min(abs(Mac.s-Mac.rw(1)));
Bs1 = B1(II,:);
Bs2 = B2(II-1,:);
Bs3 = B2(II,:);

Bs10 = spline(Mac.chi,Bs1,chis0);
Bs20 = spline(Mac.chi,Bs2,chis0);
Bs30 = spline(Mac.chi,Bs3,chis0);

chis = linspace(-pi,pi,1001);

[smin,II] = min(abs(Mac.s-Mac.rw(2)));
Bs1 = B1(II,:);
Bss1 = spline(Mac.chi,Bs1,chis);




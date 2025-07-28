% test function MacGetRatSurf.m

s = linspace(0,1,301);
q = 10*(s-0.3).^2 + 1.2;

[Ires,Sres,Qres] = MacGetRatSurf(s,q,3,[-29:29],1);


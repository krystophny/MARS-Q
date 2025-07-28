function MacResetCoil(R,Z)

global Mac

Mac.coil = Mac.coilN;
Dm = [];  km = [];
for kc = 1:size(Mac.coilN,1);
  I0 = []; J0 = []; D0 = [];
  for k=Mac.Ns1:Mac.Ns
    [Ri, Rj] = meshgrid(R(k,:),R(k,:));
    [Zi, Zj] = meshgrid(Z(k,:),Z(k,:));
    Ri = Ri*Mac.R0EXP;
    Rj = Rj*Mac.R0EXP;
    Zi = Zi*Mac.R0EXP;
    Zj = Zj*Mac.R0EXP;
    tmp = sqrt( (Ri-Mac.coilN(kc,1)).^2 + (Zi-Mac.coilN(kc,2)).^2 );
    tmp = tmp + sqrt( (Rj-Mac.coilN(kc,3)).^2 + (Zj-Mac.coilN(kc,4)).^2 );
    [Y,II] = min(tmp);
    [X,JJ] = min(Y);
    I0 = [I0; II(JJ)];
    J0 = [J0; JJ];
    D0 = [D0; X];
  end
  [Dmin,kmin] = min(D0);
  Imin = I0(kmin);
  Jmin = J0(kmin);
  Dm = [Dm; Dmin];
  km = [km; kmin];
  Mac.coil(kc,1) = Mac.s(Mac.Ns1+kmin-1);
  chi1 = Mac.chi(Jmin)/pi;
  chi2 = Mac.chi(Imin)/pi;
  if chi1 > 1.0, chi1 = chi1 - 2; end
  if chi2 > 1.0, chi2 = chi2 - 2; end
  if chi1 > chi2, tmp = chi1; chi1 = chi2; chi2 = tmp; end
  Mac.coil(kc,2) = chi1;
  Mac.coil(kc,3) = chi2;
end

CoilGeometry = [km Mac.coil Dm]



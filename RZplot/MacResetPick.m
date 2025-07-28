function MacResetPick(R,Z)

global Mac

I0 = []; J0 = []; R0=[]; Z0=[];
for kc = 1:size(Mac.pickN,1);
    Ri = R(Mac.Ns1:Mac.Ns,:)*Mac.R0EXP;
    Zi = Z(Mac.Ns1:Mac.Ns,:)*Mac.R0EXP;
    tmp = sqrt( (Ri-Mac.pickN(kc,1)).^2 + (Zi-Mac.pickN(kc,2)).^2 );
    [Y,II] = min(tmp);
    [X,JJ] = min(Y);
    I0 = [I0; II(JJ)];
    J0 = [J0; JJ];
    R0 = [R0; Ri(II(JJ),JJ)];
    Z0 = [Z0; Zi(II(JJ),JJ)];
end

PickGeometry = [I0 Mac.pickN(:,1) R0 Mac.pickN(:,2) Z0 Mac.s(Mac.Ns1+I0-1) (Mac.chi(J0))'/pi]
disp(['IPICK=' num2str(I0')])
disp(['CPICK=' num2str(Mac.chi(J0)/pi)])

if Mac.plot_pick>0
   figure(Mac.plot_pick)
   plot(Mac.pickN(:,1),Mac.pickN(:,2),'ro','LineWidth',2,'MarkerSize',9), hold on
   plot(R0,Z0,'b+','LineWidth',2,'MarkerSize',9), hold on
end


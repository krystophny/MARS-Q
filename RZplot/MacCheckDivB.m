function [DivBM] = MacCheckDivB(BM1,BM2,BM3)
%% check DivB along the surface s=rw(1), at half-integer mesh

global Mac

[Y,I0] = min(abs(Mac.s-Mac.rw(1)));
I1 = I0-1;
mm = (Mac.Mm(:))';

DivBM = (BM1(I0,:)-BM1(I1,:))/(Mac.s(I0)-Mac.s(I1));
DivBM = DivBM + i*(mm.*BM2(I1,:) + Mac.n*BM3(I1,:));

ModBM = abs((BM1(I0,:)+BM1(I1,:))/2).^2;
ModBM = sqrt(ModBM + abs(BM2(I1,:)).^2 + abs(BM3(I1,:)).^2);
ModBM = sum(ModBM)/length(ModBM);

DivBM = DivBM/ModBM;

if Mac.plot_DivB>0  

figure(Mac.plot_DivB)
plot(Mac.Mm,real(DivBM),'r-o','LineWidth',2), hold on
plot(Mac.Mm,imag(DivBM),'b--o','LineWidth',2), hold on
ylabel('DivB/|B|','FontSize',14),
xlabel('m','FontSize',14),

end




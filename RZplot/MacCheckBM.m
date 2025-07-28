function MacCheckBM(BM1,BM2,BM3,BM1C,BM2C,BM3C)

global Mac

[smin,II] = min(abs(Mac.s-Mac.rw(1)));

BM1C = rot90(BM1C);
BM2C = rot90(BM2C);
BM3C = rot90(BM3C);
BM2(II,:)
err1 = sqrt(sum(abs(BM1C-BM1(II,:)).^2)/sum(abs(BM1(II,:)).^2));
err2 = sqrt(sum(abs(BM2C-BM2(II,:)).^2)/sum(abs(BM2(II,:)).^2));
err3 = sqrt(sum(abs(BM3C-BM3(II,:)).^2)/sum(abs(BM3(II,:)).^2));

disp(['MacCheckBM: relative error=',num2str(err1),' ',num2str(err2),' ',num2str(err3),' ',num2str(err1+err2+err3)])

if Mac.plot_BMcheck>0
figure(Mac.plot_BMcheck)
subplot(3,2,1), plot(Mac.Mm,real(BM1(II,:)),'bo',Mac.Mm,real(BM1C),'r+','LineWidth',2,'MarkerSize',9), hold on,
                ylabel('Re[B_m^1]','FontSize',16),
subplot(3,2,2), plot(Mac.Mm,imag(BM1(II,:)),'bo',Mac.Mm,imag(BM1C),'r+','LineWidth',2,'MarkerSize',9), hold on,
                ylabel('Im[B_m^1]','FontSize',16),
subplot(3,2,3), plot(Mac.Mm,real(BM2(II,:)),'bo',Mac.Mm,real(BM2C),'r+','LineWidth',2,'MarkerSize',9), hold on,
                ylabel('Re[B_m^2]','FontSize',16),
subplot(3,2,4), plot(Mac.Mm,imag(BM2(II,:)),'bo',Mac.Mm,imag(BM2C),'r+','LineWidth',2,'MarkerSize',9), hold on,
                ylabel('Im[B_m^2]','FontSize',16),
subplot(3,2,5), plot(Mac.Mm,real(BM3(II,:)),'bo',Mac.Mm,real(BM3C),'r+','LineWidth',2,'MarkerSize',9), hold on,
                ylabel('Re[B_m^3]','FontSize',16),
subplot(3,2,6), plot(Mac.Mm,imag(BM3(II,:)),'bo',Mac.Mm,imag(BM3C),'r+','LineWidth',2,'MarkerSize',9), hold on,
                ylabel('Im[B_m^3]','FontSize',16),
end

function [DivJM] = MacCheckDivJ(JM1,JM2,JM3)
%% check DivJ along the surface s=rw(1), at integer mesh

global Mac

[Y,I0] = min(abs(Mac.s-Mac.rw(1)));
I0 = Mac.Iratsurf(1);

I1 = I0-1;  
mm = (Mac.Mm(:))';
h1 = (Mac.s(I0)-Mac.s(I0-1))/2;
h2 = (Mac.s(I0+1)-Mac.s(I0))/2;

DJM1 = (JM1(I0,:)-JM1(I1,:))/(h1+h2);
DivJM = DJM1 + i*(mm.*JM2(I0,:) + Mac.n*JM3(I0,:));

ModJM = abs((h1*JM1(I0,:)+h2*JM1(I1,:))/(h1+h2)).^2;
ModJM = sqrt(ModJM + abs(JM2(I0,:)).^2 + abs(JM3(I0,:)).^2);
ModJM = sum(ModJM)/length(ModJM);

a = mm;
b = DJM1;
c = JM2(I0,:);
d = JM3(I0,:);
res = [a(:) b(:) c(:) d(:)]

DivJM = DivJM/ModJM;

if Mac.plot_DivJ>0  

figure(Mac.plot_DivJ)
plot(Mac.Mm,real(DivJM),'r-o','LineWidth',2), hold on
plot(Mac.Mm,imag(DivJM),'b--o','LineWidth',2), hold on
ylabel('DivJ/|J|','FontSize',14),
xlabel('m','FontSize',14),

end




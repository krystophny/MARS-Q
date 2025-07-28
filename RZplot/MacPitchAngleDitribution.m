% plot pitch angle distribution for NBI EP models
% KNBI=13

Z0 = 0.8;        %=ZZETA0
ZS = 0.25;       %varies with magnetic surface
DZ = sqrt(0.15); %=DZETA0

C1 = [1 1 1 1 0 0]; ZI1 = [Z0 2-Z0 -Z0 -2+Z0 0 0];
C2 = [1 1 1 1 0 0]; ZI2 = [Z0 2-Z0 -Z0 -2+Z0 0 0];
C3 = [1 1 0 0 1 1]; ZI3 = [Z0 2-Z0 -Z0 -2+Z0 -2*ZS+Z0 -2*ZS+2-Z0];


Z = linspace(-1,1,101); 
f = Z*0;
for k=1:length(Z)
    if Z(k)>=ZS & Z(k)<=1,   C=C1; ZI=ZI1; end
    if Z(k)>=-ZS & Z(k)<=ZS, C=C2; ZI=ZI2; end
    if Z(k)>=-1 & Z(k)<=-ZS, C=C3; ZI=ZI3; end

    for j=1:length(C)
        f2   = exp(-(Z(k)-ZI(j))^2/DZ^2)/(2*sqrt(pi)*DZ);
        f(k) = f(k) + C(j)*f2;
    end
end

hf=figure(1);
plot(Z,f,'b-','LineWidth',2), hold on
xlabel('\zeta','FontSize',18,'FontWeight','Bold')
ylabel('f^1','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
a=axis; 
plot([-ZS -ZS],[a(3) a(4)],'k--'), 
plot([ ZS  ZS],[a(3) a(4)],'k--'), 



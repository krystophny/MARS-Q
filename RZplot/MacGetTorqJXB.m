% compute toroidal component of the jxb torque density
% at half-integer points

function [Tjxb,Ajxb] = MacGetTorqJXB(BM1,BM2,JM1,JM2,jacobian)

global Mac

II = 1:Mac.Ns1;
IX = 1:Mac.Ns1-1;

% compute q-profile at half-integer points
%q  = Mac.q;
q  = Mac.s(1:Mac.Ns1);
qx = (q(1:end-1)+q(2:end))/2;

% compute J0 for normalisation purpose
J0 = 2*pi*sum(jacobian(II,1:end-1),2)*(Mac.chi(2)-Mac.chi(1));
J0 = (J0(1:end-1)+J0(2:end))/2;
J0 = J0*ones(1,size(BM1,2));

% compute Tjxb due to Re(J1)*Re(B2)
T12R = real(JM1(IX,:)).*real(BM2(IX,:))*2*pi^2./J0;

% compute Tjxb due to Im(J1)*Im(B2)
T12I = imag(JM1(IX,:)).*imag(BM2(IX,:))*2*pi^2./J0;

% compute Tjxb due to Re(J2)*Re(B1)
T21R = real(JM2(II,:)).*real(BM1(II,:));
T21R = (T21R(1:end-1,:)+T21R(2:end,:))/2*2*pi^2./J0;

% compute Tjxb due to Im(J2)*Im(B1)
T21I = imag(JM2(II,:)).*imag(BM1(II,:));
T21I = (T21I(1:end-1,:)+T21I(2:end,:))/2*2*pi^2./J0;

% compute total Tjxb density
Tjxb = sum(T12R+T12I-T21R-T21I,2);
Tjxb(end) = Tjxb(end-1);
%Tjxb = Tjxb/max(abs(Tjxb));

s = q; 
save TorqData s J0 JM1 JM2 BM1 BM2 

% compute accumulated torque
Ajxb = zeros(Mac.Ns1,1);
Ajxb(1) = 0;
for k=1:Mac.Ns1-1
    Ajxb(k+1) = Ajxb(k) + Tjxb(k)*J0(k,1)*(Mac.s(k+1)-Mac.s(k));
end

% plot jxb torque density
if Mac.plot_Tjxb
   hf = figure(10*Mac.plot_Tjxb+1);
   MN = Mac.Mm;
   MM = MN - Mac.Mm(1) + 1;
   Y  = T12R(:,MM); Cn=max(max(abs(Y))); %Y = Y/Cn;
   hp = plot(qx,Y,'-','LineWidth',1); hold on,
   ylabel('T_{jxb}^{12R}','FontSize',18,'FontWeight','Bold')
   xlabel('q','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(qx(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end

   hf = figure(10*Mac.plot_Tjxb+2);
   MN = Mac.Mm;
   MM = MN - Mac.Mm(1) + 1;
   Y  = T12I(:,MM); Cn=max(max(abs(Y))); %Y = Y/Cn;
   hp = plot(qx,Y,'-','LineWidth',1); hold on,
   ylabel('T_{jxb}^{12I}','FontSize',18,'FontWeight','Bold')
   xlabel('q','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(qx(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end

   hf = figure(10*Mac.plot_Tjxb+3);
   MN = Mac.Mm;
   MM = MN - Mac.Mm(1) + 1;
   Y  = T21R(:,MM); Cn=max(max(abs(Y))); %Y = Y/Cn;
   hp = plot(qx,Y,'-','LineWidth',1); hold on,
   ylabel('T_{jxb}^{21R}','FontSize',18,'FontWeight','Bold')
   xlabel('q','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(qx(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end

   hf = figure(10*Mac.plot_Tjxb+4);
   MN = Mac.Mm;
   MM = MN - Mac.Mm(1) + 1;
   Y  = T21I(:,MM); Cn=max(max(abs(Y))); %Y = Y/Cn;
   hp = plot(qx,Y,'-','LineWidth',1); hold on,
   ylabel('T_{jxb}^{21I}','FontSize',18,'FontWeight','Bold')
   xlabel('q','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(qx(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end

   hf = figure(10*Mac.plot_Tjxb+5);
   Tjxb = Tjxb*Mac.B0EXP^2/4e-7/pi;
   sx   = qx;
   hp = plot(sx,Tjxb,'-','LineWidth',1); hold on,
   save DenJXB sx Tjxb
   ylabel('T_{jxb}','FontSize',18,'FontWeight','Bold')
   xlabel('q','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

   hf = figure(10*Mac.plot_Tjxb+6);
   hp = plot(q,Ajxb,'k-','LineWidth',2); hold on,
   ylabel('\int{T}_{jxb}dr','FontSize',18,'FontWeight','Bold')
   xlabel('q','FontSize',18,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   %adds patch
   if 1==0
   axis([1.5 2.65 -0.06 0.01])
   xx = [1.5 1.9534 1.9534     2.1019     2.1019 2.2662 2.2662  2.4     2.4     2.65];
   yy = [0.0 0.0   -2.0167e-3 -2.0167e-3 -0.031 -0.031 -0.0416 -0.0416 -0.0499 -0.0499];
   plot(xx,yy,'k--','LineWidth',1)
   end
end

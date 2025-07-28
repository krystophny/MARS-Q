function [JparaM2,JparaM3]=MacTestStarEqLi(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian,TEQ,DPSIEQ,J1,J2,J3)

global Mac  SDIR

N = Mac.Ns1;

%============ RS: Part II =================================================
G12   = dRds(1:N,:).*dRdchi(1:N,:) + dZds(1:N,:).*dZdchi(1:N,:);
G22   = dRdchi(1:N,:).^2 + dZdchi(1:N,:).^2;
TEQN  = TEQ*ones(1,Mac.Nchi);  
DPSN  = DPSIEQ*ones(1,Mac.Nchi);
BEQ   = sqrt(DPSN.^2.*G22./jacobian(1:N,:).^2+TEQN.^2./R(1:N,:).^2);
BEQ(1,:) = BEQ(2,:);
%JparaN2  = ( (G12.*J1(1:N,:)+G22.*J2(1:N,:)).*DPSN./jacobian(1:N,:))./BEQ; 
%JparaN3  = ( TEQN.*J3(1:N,:) )./BEQ; 
JparaN2  = (G12.*J1(1:N,:)+G22.*J2(1:N,:))./jacobian(1:N,:); 
JparaN3  = R(1:N,:).^2.*J3(1:N,:)./jacobian(1:N,:); 

expmchi = exp(-Mac.chi'*Mac.Mm'*i);
JparaM2 = JparaN2*expmchi*(Mac.chi(2)-Mac.chi(1))/2/pi;
JparaM3 = JparaN3*expmchi*(Mac.chi(2)-Mac.chi(1))/2/pi;

if 1==0
%save Li JparaM2 -ascii

hf=figure(10*Mac.plot_testStarEq+0);
%MN = Mac.Mm;
MN = 2;
MM = MN - Mac.Mm(1) + 1;

   ax = [0.9 0.906];
Y = real(JparaM2(:,MM)); Cn=max(max(abs(Y))); %Y = Y/Cn;
hp=plot(Mac.s(1:N),Y,'-','LineWidth',1); hold on,
ylabel('Re[J_2]','FontSize',18,'FontWeight','Bold')
%xlabel('s\equiv\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
xlabel('s','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes');
set(ha,'FontSize',16,'FontWeight','Bold')
for k=1:length(MN)
    c = get(hp(k),'Color');
    [X,I]=max(abs(Y(:,k)));
    text(Mac.q(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
end
a= axis;
for k=1:length(Mac.Iratsurf)
   plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
end
 xlim(ax);

hf=figure(10*Mac.plot_testStarEq+1);
%MN = Mac.Mm;
MN = 2;
MM = MN - Mac.Mm(1) + 1;
Y = imag(JparaM3(:,MM)); %Cn=max(max(abs(Y))); Y = Y/Cn;
hp=plot(Mac.s(1:N),Y,'-','LineWidth',1); hold on,
ylabel('Im[J_3]','FontSize',18,'FontWeight','Bold')
xlabel('s\equiv\psi_p^{1/2}','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes');
set(ha,'FontSize',16,'FontWeight','Bold')
for k=1:length(MN)
    c = get(hp(k),'Color');
    [X,I]=max(abs(Y(:,k)));
    text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
end
a= axis;
for k=1:length(Mac.Iratsurf)
    plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
end
 xlim(ax);

% ==== plot RM ===
   hf=figure(10*Mac.plot_testStarEq+2);
 
   MN = 2;
   MM = MN - Mac.Mm(1) + 1;
   Rm = JparaM2(:,MM) - (MN/Mac.n) * JparaM3(:,MM);
   ax = [0.9 0.906];

   subplot(2,2,1);
   Y = real(Rm); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(1:N),Y,'LineWidth',2); hold on,
   ylabel('real(Rm)','FontSize',16,'FontWeight','Bold') 
   xlabel('s\equiv\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   xlim(ax);

   subplot(2,2,2);
   Y = imag(Rm); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(1:N),Y,'LineWidth',2); hold on,
   ylabel('imag(Rm)','FontSize',18,'FontWeight','Bold') 
   xlabel('s\equiv\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   xlim(ax);
  
   subplot(2,2,3);
   Y = abs(Rm); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(1:N),Y,'LineWidth',2); hold on,
   ylabel('abs(Rm)','FontSize',16,'FontWeight','Bold') 
   xlabel('s\equiv\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   xlim(ax);
   
 % ==== plot RSII = eta* RM ===
   hf=figure(10*Mac.plot_testStarEq+3);

   II = 1:Mac.Ns1;
   MN = 2;
   MM = MN - Mac.Mm(1) + 1;
   
   %MacInputPROFEQ('PROFEQ.OUT');
   
   RSII = (Mac.Resieta(II)).*Rm(II);
   %RSII = 1.0e-9*Rm(II);
   %ax = [0.9 0.906];
%    ax = [0.82 0.95];
      ax = [0.9 0.906];
   
%    YY(:,1) = real(Rm);
%    YY(:,2) = imag(Rm);
%    YY(:,3) = real(RSII);
%    YY(:,4) = imag(RSII);
%    YY(:,5) = Mac.Resieta;
% save Li1 YY -ascii


   subplot(2,2,1);
   Y = real(RSII); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(1:N),Y,'LineWidth',2); hold on,
   ylabel('real(RSII)','FontSize',16,'FontWeight','Bold') 
   xlabel('s\equiv\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   xlim(ax);

   subplot(2,2,2);
   Y = imag(RSII); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(1:N),Y,'LineWidth',2); hold on,
   ylabel('imag(RSII)','FontSize',18,'FontWeight','Bold') 
   xlabel('s\equiv\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   xlim(ax);
  
   subplot(2,2,3);
   Y = abs(RSII); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(1:N),Y,'LineWidth',2); hold on,
   ylabel('abs(RSII)','FontSize',16,'FontWeight','Bold') 
   xlabel('s\equiv\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   title('RSII = \eta * Rm');
   xlim(ax);
%========================================================================== 
% ====== RS:Part: I ===================
   XPLASMA = load([SDIR 'XPLASMA.OUT']);
   Mac.Nm1 = XPLASMA(1,1);
   Mac.n  = round(XPLASMA(1,3));
   Mac.Mm = round(XPLASMA(2:Mac.Nm1+1,1));
   VM1 = XPLASMA(Mac.Nm1+2+Mac.Ns1:end,1) + XPLASMA(Mac.Nm1+2+Mac.Ns1:end,2)*i;
   %VM1 = reshape(VM1,Mac.Ns1,Mac.Nm1)*Mac.VNORM;
   VM1 = reshape(VM1,Mac.Ns1,Mac.Nm1);

% ==== PLOT =================
   hf=figure(10*Mac.plot_testStarEq+4);
   %MN = Mac.Mm;
   II = 1:Mac.Ns1;
   MN = 2;  
   MM = MN - Mac.Mm(1) + 1;
%    ax = [0.82 0.95];
      ax = [0.9 0.906];
   
   Y1 = Mac.dpsi(II).*(MN + Mac.n*Mac.q(II)).* Mac.omega(II).* VM1(:,MM);
   subplot(2,2,1);
   Y = real(Y1);
   %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(II),Y,'-','LineWidth',2); hold on,
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Re[RSI]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   xlim(ax);
   
   subplot(2,2,2);
   Y = imag(Y1);
   %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(II),Y,'-','LineWidth',2); hold on,
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Im[RSI]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   xlim(ax);
   
   subplot(2,2,3);
   Y = abs(Y1); 
   %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(II),Y,'-','LineWidth',2); hold on,
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('|RSI|','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   title('RSI = dpsi*(m-nq)*omega*XM1');
   xlim(ax);

%==========================================================================   
    % ==== plot RS = RSI + RSII ===
   hf=figure(10*Mac.plot_testStarEq+5);

%    ax = [0.82 0.95];
      ax = [0.9 0.906];
   Y2 = Y1 + RSII;

   subplot(2,2,1);
   Y = real(Y2); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(1:N),Y,'LineWidth',2); hold on,
   ylabel('real(RS)','FontSize',16,'FontWeight','Bold') 
   xlabel('s\equiv\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   xlim(ax);

   subplot(2,2,2);
   Y = imag(Y2); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(1:N),Y,'LineWidth',2); hold on,
   ylabel('imag(RS)','FontSize',18,'FontWeight','Bold') 
   xlabel('s\equiv\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   xlim(ax);
  
   subplot(2,2,3);
   Y = abs(Y2); %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(1:N),Y,'LineWidth',2); hold on,
   ylabel('abs(RS)','FontSize',16,'FontWeight','Bold') 
   xlabel('s\equiv\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   title('RS = RSI + RSII');
   xlim(ax);
   
%==========================================================================
% ====== LS:Part ======================
   BPLASMA = load([SDIR 'BPLASMA.OUT']);
   Mac.Nm1 = BPLASMA(1,1);
   Ns = BPLASMA(1,2);
   BM1 = BPLASMA(Mac.Nm1+2:end,1) + BPLASMA(Mac.Nm1+2:end,2)*i;
   BM1 = reshape(BM1,Ns,Mac.Nm1);
   BM1 = BM1(1:Mac.Ns,:);
   %BM1 = BM1*Mac.BNORM;  %[T]

% ==== PLOT ======
   hf=figure(10*Mac.plot_testStarEq+6);

   %MN = Mac.Mm;
   II = 1:Mac.Ns1;
   MN = 2;  
   MM = MN - Mac.Mm(1) + 1;
   
   Y1 = Mac.omega(II).* BM1(II,MM);
%    ax = [0.82 0.95];
      ax = [0.9 0.906];
   subplot(2,2,1);
   Y = real(Y1);
   %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(II),Y,'r-','LineWidth',2); hold on,
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Re[LS]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   xlim(ax);
   %ylim([-1e-8 2e-8]);
   
   subplot(2,2,2);
   Y = imag(Y1);
   %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(II),Y,'r-','LineWidth',2); hold on,
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Im[LS]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   xlim(ax);
   %ylim([-1e-8 2e-8]);
   
   subplot(2,2,3);
   Y = abs(Y1); 
   %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(II),Y,'r-','LineWidth',2); hold on,
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('|LS|','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   title('LS = omega*BM1');
   xlim(ax);
   %ylim([-1e-8 2e-8]);

   % ==== PLOT ======
   hf=figure(10*Mac.plot_testStarEq+7);

   %MN = Mac.Mm;
   II = 1:Mac.Ns1;
   MN = 2;  
   MM = MN - Mac.Mm(1) + 1;
   
   Y1 = BM1(II,MM);
%    ax = [0.902 0.904];
      ax = [0.9 0.906];
   subplot(2,2,1);
   Y = real(Y1);
   %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(II),Y,'-','LineWidth',2); hold on,
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Re[BM1]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   xlim(ax);
   %axis([0.90 0.906 -1.0e-3 1.0e-3]);
   
   subplot(2,2,2);
   Y = imag(Y1);
   %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(II),Y,'-','LineWidth',2); hold on,
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Im[BM1]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   xlim(ax);
   
   subplot(2,2,3);
   Y = abs(Y1); 
   %Cn=max(max(abs(Y))); Y = Y/Cn;
   hp=plot(Mac.s(II),Y,'-','LineWidth',2); hold on,
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('|BM1|','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(MN)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(MN(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   a= axis;
   for k=1:length(Mac.Iratsurf)
       plot([Mac.s(Mac.Iratsurf(k)) Mac.s(Mac.Iratsurf(k))],[a(3) a(4)],'k--')
   end
   plot(ax,[0 0],'k--');
   title('LS = omega*BM1');
   xlim(ax);
   
end

end

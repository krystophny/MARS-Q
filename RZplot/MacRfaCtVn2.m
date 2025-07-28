%function MacRfaCtBn2
%% The 2nd (fast) method to
%% transform from X1M in MARS coordinate to X1M in PEST coordinate
%% and plot the BnM poloidal harmonics
%% before running this subroutine, make sure that 
%% MacMain* is a procedure, not a function
%%   1) RMZM_F & XPLASMA are used in the MacMain* subroutine
%%   2) SDIR defined here agrees with those in the MacMain* subroutine 
%%   3) files RMZM_F_EQAC, RMZM_F_PEST, XPLASMA(EQAC) exist in SDIR
%%   4) comment out Mac.RunB in MacMain*
%%   5) FEEDI is explicitely assumed to be =1 in MacMain*
%% NB: BnPEST is actually B1PEST=J*(b\cdot\nabla s). Only B1PEST has the property 
%%     of vanishing amplitude at the rational surfaces in ideal MHD

global Mac 
global SDIR 

%SDIR = '/home/liuy/Work/MAST-U45272/Data/'; facn=1.0; SMAIN='D3D_ALL_B';
%SDIR = '/home/liuy/Work/DIII-D/D3D_QH_CETOP/Data/163520.1750/n2/'; facn=1.0; SMAIN='D3D_ALL_B';
SDIR = '/home/liuy/Work/NSTX/NSTX132543/Data/'; facn=1.0; SMAIN='D3D_ALL_B';
  
facn      = 1.0;
Mac.RunJ  = 0;
Mac.RunB  = 0;
n         = 19;
CheckQ    = 0;
CheckErgos= 0;
PlotBn    = 6;
PlotB1E   = 0;
PlotQ     = 5;
PlotPhi   = 0;
LSS       = '-o'; 
%LSC       = [(kclr-1)/(Nclr-1) (kclr-1)*(Nclr-kclr)*4/Nclr^2 (Nclr-kclr)/(Nclr-1)]
LSC       = [0 1 0];

mk = -9:170;
kxp= 2;
mi = mk - mk(1) + 1;
Mac.edge = 0.999;
Mac.core = 0.21;
sedge = Mac.edge;
  
eval(['copyfile ' SDIR 'RMZM_F_EQAC ' SDIR 'RMZM_F.OUT']);
Mac.RunV = 1;
eval(['MacMain' SMAIN]);
II = 1:Mac.Ns1;
BnEQAC = Vn(II,:);
R_EQAC = R(II,:);
Z_EQAC = Z(II,:);

%compute Bn vs. geometrical angle theta at the plasma surface
IED = find(Mac.s<sedge);
IED = IED(end);
Rs = R_EQAC(IED,:);
Zs = Z_EQAC(IED,:);
Rc = (min(Rs)+max(Rs))/2;
Zc = (min(Zs)+max(Zs))/2;
Tg = atan2(Zs-Zc,Rs-Rc);
[Y,IXU] = max(Zs); 
[Y,IXL] = min(Zs); 
[Y,IXM] = max(Rs); 
BnEDGE = BnEQAC(IED,:);
TgX    = Tg(IXL);

eval(['copyfile ' SDIR 'RMZM_F_PEST ' SDIR 'RMZM_F.OUT']);
Mac.RunV = 0;
eval(['MacMain' SMAIN]);
ss     = Mac.s(II);
R_PEST = R(II,:);
Z_PEST = Z(II,:);
G22_PEST  = dRdchi(II,:).^2 + dZdchi(II,:).^2;  G22_PEST(1,:) = G22_PEST(2,:);

BnPEST  = griddata(R_EQAC,Z_EQAC,BnEQAC,R_PEST,Z_PEST);
BnPEST = BnPEST.*sqrt(G22_PEST).*R_PEST./jacobian(II,:);  %X1PEST

%save and plot Br along PEST poloidal angle at one s-surface
if 1==0
  [sp,Ip]=min(abs(Mac.s-0.9))
  res = [Mac.chi' real(BnPEST(Ip,:))'*Mac.B0EXP*4e+3/Mac.R0EXP  imag(BnPEST(Ip,:))'*Mac.B0EXP*4e+3/Mac.R0EXP];
  save B1.txt res -ascii
  figure(90)
  plot(res(:,1),res(:,2),'r-',res(:,1),res(:,3),'b--')  
end

ss = Mac.s(1:Mac.Ns1);
BnPEST = BnPEST(1:Mac.Ns1,:);


expmchi = exp(-Mac.chi'*mk*i);
BMnPEST = BnPEST*expmchi*(Mac.chi(2)-Mac.chi(1))/2/pi;

%patch possible NaN's of BMnPEST near the plasma boundary
if 1==0
for k=1:size(BMnPEST,2)
    INOR = find(isnan(BMnPEST(:,k))==0);
    if length(INOR) < size(BMnPEST,1)
       BMnPEST(:,k) = pchip(ss(INOR),BMnPEST(INOR,k),ss);
    end
end
end

%patch BnEQAC near plasma edge
if 1==1
II = find(Mac.s<=sedge);
II = II(end);
for k=II+1:Mac.Ns1
    BnEQAC(k,:) = BnEQAC(II,:);
end
end

%modify phase of BMnPEST
if 1==0
m2 = 2-mk(1)+1;
[Y,I2]=max(abs(imag(BMnPEST(:,m2))));
a = BMnPEST(I2,m2);
p = phase(a);
BMnPEST = BMnPEST*exp(-i*p);
end

if 1==1
  save BnMat ss mk Tg BMnPEST BnEDGE
  eval(['movefile ' 'BnMat.mat ' SDIR]);

  B1E = BnPEST(end,:); B1E=B1E(:);
  res = [Mac.chi(:) real(B1E) imag(B1E)];
  save FootPrintData_B1E res -ascii  
else
  eval(['load ' SDIR 'BnMat.mat']);
end

mm = mk;
II = mm - mk(1) + 1;
mk = mk(II); 
BnPEST2 = BnPEST;
BnPEST  = BMnPEST(:,II)/facn;
BnEDGE  = BnEDGE/facn;
BnPEST(1,:) = BnPEST(2,:);

if CheckQ > 0
   dataq = load([SDIR 'PROFEQ.OUT']);
   s = dataq(:,1);
   q = dataq(:,2);
   dpsi = dataq(:,12);

   % normalisation according to ERGOS
   if CheckErgos>0
      if CheckErgos==1, fac_ERGOS = q/4.*dpsi*Mac.R0EXP^2*Mac.B0EXP*1e-4; end
      if CheckErgos==2, fac_ERGOS = q.*dpsi*Mac.R0EXP; end
      if CheckErgos==3, fac_ERGOS = q/Mac.B0EXP; end
      for k=1:size(BnPEST,2)
          BnPEST(2:end,k) = BnPEST(2:end,k)./fac_ERGOS(2:end);
      end
   end
end


if PlotBn > 0
   if 1==0
   hf=figure(10*PlotBn + 0);
   res = abs(BnPEST);
   II = find(ss(:,1)<1.5);
   pcolor(mk,ss(II,:).^2,res(II,:)), hold on, shading interp
   %axis([min(mk) max(mk) 0 1])
   axis([-20 20 0 1])
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('\psi_p','FontSize',16,'FontWeight','Bold')
   title('|\xi^1|','FontSize',16,'FontWeight','Bold')
   colorbar,  colormap(hot)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
   %map = colormap; N = size(map,1); a = linspace(0,1,N)';  b = a.^0.5*ones(1,3); colormap(map.*b)
    end

   %plot rational surfaces
   mm3 = Mac.mm_plot; II3=mm3-mk(1)+1;
   if CheckQ > 0
      mq = [ceil(min(q)*abs(n)):max(mk)];
      qq = mq/abs(n);
      [sq,qn] = MacFindX(s,q,qq);
      mq = qn*n;
      %plot(mq,sq,'b+','LineWidth',3,'MarkerSize',8), hold on,
      save Qsurf mq sq
      mm2 = mq; II2=mm2-mk(1)+1;
      mm3 = [mq(end)+1:mq(end)+6]; II3=mm3-mk(1)+1; 
     %mm3 = mq(2); II3=mm3-mk(1)+1;
     %mm3 = mk; II3=mm3-mk(1)+1;
     %II3 = find(mk==11);  mm3=mk(II3);
     %II3 = II2;
   end

   if CheckQ > 0 & 1==0
      hf=figure(10*PlotQ + 0);
      plot(s.^2,q,LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
      xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
      ylabel('q','FontSize',16,'FontWeight','Bold')
      ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
      a = axis;
      for k=1:length(sq)
          plot([0 1],[mq(k) mq(k)]/n,'k--'), hold on,
          plot([sq(k)^2 sq(k)^2],[a(3) a(4)],'k--'), hold on,
      end
   end
   
   %plot Bn near the plasma boundary
   if 1==0
   hf=figure(10*PlotBn + 1);
   s2 = ss.^2;
   s0 = 0.9;
   JJ = find(s2>s0);
   pcolor(mk,s2(JJ),abs(BnPEST(JJ,:))), hold on, shading interp
   contour(mk,s2(JJ),abs(BnPEST(JJ,:)),'k-'), hold on,
   axis([min(mk) max(mk) s0 1])
	xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('\psi_p','FontSize',16,'FontWeight','Bold')
   title('|\xi^1| [mm]','FontSize',16,'FontWeight','Bold')
   colorbar,  colormap(hot)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
   if CheckQ > 0
      k = find(sq>s0);
      plot(mq(k),sq(k).^2,'b+','LineWidth',3,'MarkerSize',8), hold on,
   end
   end
    
   if 1==0
   hf=figure(10*PlotBn + 2);
   plot(ss.^2,real(BnPEST(:,II3)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('Re(\xi^1_m) [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq].^2,[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 
   end

   if 1==0
   hf=figure(10*PlotBn + 3);
   plot(ss.^2,imag(BnPEST(:,II3)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('Im(\xi^1_m) [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq].^2,[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 
   end

   hf=figure(10*PlotBn + 4);
   II  = find(Mac.s<Mac.core);
   II  = II(end):Mac.Ns1;
   II4 = II3;
   %II4 = [II3(1) II3(4)];
   [D,I] = max(abs(BnPEST(II,II4)),[],1);
   [C,J] = max(D);
   Pm1 = angle(BnPEST(I(J),II4(J)));
   BnPEST = BnPEST*exp(-i*Pm1);
   %C   = 1;
   for ks=1:length(II4)
   %hs=subplot(2,2,ks,'align');
   Y=real(BnPEST(:,II4(ks)))/C;
   hp=plot(ss(II).^kxp,Y(II),LSS(1:end-1),'LineWidth',2); hold on,
   xlabel('s=\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   if kxp==2
      xlabel('\psi_p=s^2','FontSize',16,'FontWeight','Bold')
   end
   ylabel('\xi\cdot\nabla{s} [a.u.]','FontSize',16,'FontWeight','Bold')
   %ylabel('\xi\cdot\nabla{s} [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   %set(hs,'FontSize',16,'FontWeight','Bold')

   c = get(hp,'Color');
   [X,I]=max(Y);
   %text(Mac.s(I).^kxp,Y(I),int2str(II4(ks)+mk(1)-1),'FontSize',18,'FontWeight','Bold','Color',c)

   IKINK = find(ss.^2<0.5);
   IPEEL = find(ss.^2>0.8&ss.^2<sedge^2);
   FIN_RES_Akink = max(max(abs(BnPEST(IKINK,:))))*1e+3;
   FIN_RES_Apeel = max(max(abs(BnPEST(IPEEL,:))))*1e+3;
   %IKINK = find(mk<5&mk>0);
   %IPEEL = find(mk>10);
   %FIN_RES_Akink = max(max(abs(BnPEST(:,IKINK))))*1e+3;
   %FIN_RES_Apeel = max(max(abs(BnPEST(:,IPEEL))))*1e+3;
   end
   a=axis; axis([ss(II(1)).^kxp 1 a(3) a(4)])

   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      for ks=1:length(mm2)
      sqn = sq(round(mm2(ks)-mq(1)+1));
      plot([sqn; sqn].^kxp,[a(3)*ones(size(sqn)); a(4)*ones(size(sqn))],'k--'), hold on,
      end   
   end 
   %saveas(hf,'x1m1d.fig')

      %find abs(BnPEST) at rational surfaces
      if CheckQ > 0
      res = [sq; sq]; resv = sq; peaka=sq;
      for k=1:length(sq)
          [X,ISQ] = min(abs(ss-sq(k)));
          res(:,k)  = [real(BnPEST(ISQ,mq(k)-mk(1)+1)); imag(BnPEST(ISQ,mq(k)-mk(1)+1))];
          peaka(k)= max(abs(BnPEST(:,mq(k)-mk(1)+1)));
      end
      BnPEST_RS = [mq; sq; peaka; res]'
      dpsiq     = spline(s,dpsi,sq);




      if 1==0
      hf=figure(10*PlotBn + 11);
      plot(sq,abs(res(1,:)+res(2,:)*i),LSS,'LineWidth',2,'Color',LSC,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
      xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
      ylabel('|\xi^1_{res}| [mm]','FontSize',16,'FontWeight','Bold')
      ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
      end
      end

%  adjust q profile
   if CheckQ > 0 & 1==0
      sqn = sq;
      III = [-9:9];
      for k=1:length(sq)
          [X,ISQ] = min(abs(ss-sq(k)));
          [X,IQN] = min(abs(BnPEST(ISQ+III,mq(k)-mk(1)+1)));
          sqn(k)  = ss(ISQ+III(IQN));
      end
      figure(1)
      plot(s,q,'b-',sq,mq/abs(n),'b+'), hold on,
      plot(sqn,mq/abs(n),'r+'), hold on,
      II1 = find(ss<0.95);
      II2 = find(sqn>0.95); 
      sn = [ss(II1); sqn(II2)'; 1];
      qn = [q(II1);  mq(II2)'/abs(n); 9.8];
      qn = spline(sn,qn,ss);
      plot(ss,qn,'r-'), hold on,
      res = [ss qn];
      save PROFQ_EQAC res -ascii -double
      eval(['!mv ' 'PROFQ_PEST ' SDIR 'PROFQ_PEST']);
  end 

   if CheckQ > 0 & 1==0
   hf=figure(10*PlotBn + 5);
   JJ = mq - mk(1) + 1;
   plot(mk(JJ),1e+3*max(abs(BnPEST(:,JJ)),[],1),LSS,'Color',LSC,'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
   plot(mk,1e+3*max(abs(BnPEST),[],1),LSS,'Color',LSC,'LineWidth',2,'MarkerSize',7), hold on,
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('max|\xi^1_m(\psi_p)| [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   end

%  total Bn along plasma surface
   if 1==0
   hf=figure(10*PlotBn + 7);
   [x,I]=sort(Tg); 
   plot(Tg(I)*180/pi,1e+3*real(BnEDGE(I)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('geometric \theta','FontSize',16,'FontWeight','Bold')
   ylabel('Re[\xi_n] [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   end
 
%  total Xn along plasma surface
   hf=figure(10*PlotBn + 8);
   [x,I]=sort(Tg); 
   plot(Tg(I)*180/pi,abs(BnEDGE(I)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('geometric \theta','FontSize',16,'FontWeight','Bold')
   ylabel('|\xi_n| [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   a=axis;
   %plot([TgX TgX]*180/pi,[a(3) a(4)],'g-','LineWidth',6), hold on,

   res1 = Tg(I)*180/pi;
   res2 = 1e+3*abs(BnEDGE(I));
   res = [res1(:) res2(:)];
   save XnSurfGeom.txt res -ascii -double

   X = Tg(I)*180/pi;  Y = abs(BnEDGE(I));
   [Ymin,IKINK] = min(abs(X-0));
 %IPEEL = find(X<X(Ixpt));
   IPEEL = find(X<-45);
 %IPEEL = find(abs(X)>45);
   FIN_RES_DX = max(Y(IPEEL));
   FIN_RES_DM = Y(IKINK);

   if 1==0
   %plot Xn directly from MARS-F output
   datax = load([SDIR 'PROFDISP.OUT']);
   tt    = datax(:,1)*180/pi;
   XnM   = datax(:,2)+datax(:,3)*i;
   II    = find(tt>180);
   tt(II)= tt(II)-360;
   [tt,II]=sort(tt);
   XnM    = XnM(II);
   plot(tt,abs(XnM)*Mac.R0EXP*1e+3,'k-'), hold on,
   end

   % field spectrum near the plasma surface
   if CheckQ > 0 & PlotB1E>0
   IN = isnan(abs(BnPEST(:,1))); 
   II = find(IN==0);
   II = find(ss<=0.995);
   IE = II(end);
   JJ = mq - mk(1) + 1;
   hf=figure(10*PlotB1E + 1);
   plot(mk,1e+3*real(BnPEST(IE,:)),LSS,'LineWidth',1,'MarkerSize',7,'Color',LSC), hold on,
   plot(mk(JJ),1e+3*real(BnPEST(IE,JJ)),LSS,'Color',LSC,'LineWidth',1,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('|\xi^1_m(\psi_p=1)| [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   end

   %  plot B1 at plasma surface, along PEST poloidal angle 
   if PlotB1E>0
   hf=figure(10*PlotB1E + 2);
   plot(Mac.chi,1e+3*real(B1E),'-',Mac.chi,imag(B1E),'--','LineWidth',2,'Color',LSC), hold on,
   xlabel('\theta_{PEST}','FontSize',16,'FontWeight','Bold')
   ylabel('\xi^1 [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   legend('real','imaginary')
   end


%  save edge displacement on (R,Z)-grid
   res = [R_EQAC(Mac.Ns1,:)'*Mac.R0EXP Z_EQAC(Mac.Ns1,:)'*Mac.R0EXP real(BnEDGE(:))*1e+3 imag(BnEDGE(:))*1e+3];
   save EdgeDisp.txt res -ascii -double 

%plot X1 in (R_PEST,Z_PEST) space
   if 1==0
   hf=figure(10*PlotBn + 12);
   res = real(BnPEST2);
   pcolor(R_PEST*Mac.R0EXP,Z_PEST*Mac.R0EXP,res), hold on, shading interp
   xlabel('R [m]','FontSize',16,'FontWeight','Bold')
   ylabel('Z [m]','FontSize',16,'FontWeight','Bold')
   title('|\xi^1|','FontSize',16,'FontWeight','Bold')
   colorbar,  colormap(hot)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
   axis equal
   end

%plot Xn in (R,Z) space
   hf=figure(10*PlotBn + 13);
   res = abs(BnEQAC);
   II  = find(Mac.s<Mac.core);
   II = II(end):Mac.Ns1;
   res = res(II,:);
   res = res/max(max(res));
   res = log10(res);
   pcolor(R_EQAC(II,:)*Mac.R0EXP,Z_EQAC(II,:)*Mac.R0EXP,res), hold on, shading interp
   colorbar,  colormap(jet)
   axis equal
   %axis([1 2.5 -1.0 1.0])
   axis([1 2.5 -1.2 1.2])
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
   title('log10(|\xi_n|) [a.u.]','FontSize',14,'FontWeight','Bold')
   %title('log_{10}(|\xi_n|) [mm]','FontSize',14,'FontWeight','Bold')
   ylabel('Z [m]','FontSize',16,'FontWeight','Bold')
   xlabel('R [m]','FontSize',16,'FontWeight','Bold')
   %saveas(hf,'xn2d.fig')
end

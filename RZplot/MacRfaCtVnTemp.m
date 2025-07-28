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

%SDIR = '/.automount/funsrv1/root/home/yliu/D3D_MATT/';        facn=pi/2;        SMAIN='D3D134234'; %n=3
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST020333/New/';  facn=pi/2*1.7476; SMAIN='MAST020333';
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST_ALL/ResFootPrint/';  facn=pi/2*1.7476; SMAIN='MAST_ALL';
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST_ADD/';  facn=pi/2*1.7476; SMAIN='MAST_ALL';  %n=3 with 6 coils
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST_ADD/';  facn=1.9413; SMAIN='MAST_ALL';  %n=3 with 12 coils
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST_ADD/';  facn=2.1642; SMAIN='MAST_ALL';  %n=4 with 12 coils
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST_ADD/';  facn=1.6737; SMAIN='MAST_ALL';  %n=6 with 12 coils
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST_ADD/';  facn=1.0000; SMAIN='MAST_ALL';  %n=3 from MARSQ runs
%SDIR = '/.automount/funsrv4.ccfe.ac.uk/root/home1/yliu/MAST_ADD/';  facn=1.0; SMAIN='MAST_ALL_B'; 
%SDIR = '/home/yliu/MAST_ADD/';  facn=1.0; SMAIN='MAST_ALL_B'; 
SDIR = '/home/yliu/Temp/';  facn=1.0; SMAIN='AUG_B'; 
%SDIR = '/home/yliu/ITER/Scen2_V03/ResDATA_v1.2/n=3B/COMMON/';  facn=1.0; SMAIN='_ITER_B';
%SDIR = '/home/yliu/ITER/Scen2_V04/Work/';  facn=1.0; SMAIN='ITER_7KWMBX'; 
%SDIR = '/.automount/funsrv4.ccfe.ac.uk/root/home1/yliu/Scen2_V03/';  facn=1.0; SMAIN='ITER_ABT4ZL'; 
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST_ADD/';  facn=1.0; SMAIN='MAST_ALL';  %n=0
%SDIR = '/.automount/funsrv1/root/home/yliu/ANATOR/';  facn=1.0; SMAIN='ANATOR';
%SDIR = '/.automount/funsrv1/root/home/yliu/Rfx/ResEQ4b/';  facn=1.0; SMAIN='Rfx';
%SDIR = '/.automount/funsrv1/root/home/yliu/D3D117327/';       facn=pi/2;        SMAIN='D3D117327';  %n=3
%SDIR = '/.automount/funsrv1/root/home/yliu/D3D138593/Work/';       facn=1.0;        SMAIN='D3D138593';  %n=3
%SDIR = '/.automount/funsrv1/root/home/yliu/Scen2_V01/';  facn=1.0; SMAIN='ITER_15MA';
%SDIR = '/.automount/funsrv1/root/home/yliu/D3D139571/Rfa/';  facn=pi/3; SMAIN='D3D139571'; %n=1
%SDIR = '/.automount/funsrv1/root/home/yliu/JET77329/';  facn=1.0; SMAIN='JET77329';
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST22422/';  facn=1.0; SMAIN='MAST22422';
%SDIR = '/home/yliu/D3D126006/';  facn=1.0; SMAIN='D3D126006'; 

facn      = 1.0;
Mac.RunJ  = 0;
Mac.RunB  = 0;
n         = 1;
CheckQ    = 1;
CheckErgos= 0;
PlotBn    = 6;
PlotB1E   = 0;
PlotQ     = 5;
PlotPhi   = 0;
LSS       = '--d'; 
%LSC       = [(kclr-1)/(Nclr-1) (kclr-1)*(Nclr-kclr)*4/Nclr^2 (Nclr-kclr)/(Nclr-1)]
LSC       = [1 0 0];

mk = -29:29;
mi = mk - mk(1) + 1;
  
eval(['!cp ' SDIR 'rmzm_geom ' SDIR 'RMZM_F.OUT']);
Mac.RunV = 1;
eval(['MacMain' SMAIN]);
II = 1:Mac.Ns1;
BnEQAC = Vn(II,:);
R_EQAC = R(II,:);
Z_EQAC = Z(II,:);

%compute Bn vs. geometrical angle theta at the plasma surface
IED = find(Mac.s<0.995);
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

eval(['!cp ' SDIR 'rmzm_pest ' SDIR 'RMZM_F.OUT']);
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

if 1==1
  save BnMat ss mk Tg BMnPEST BnEDGE
  eval(['!mv ' 'BnMat.mat ' SDIR]);

  B1E = BnPEST(end,:); B1E=B1E(:);
  res = [Mac.chi(:) real(B1E) imag(B1E)];
  save FootPrintData_B1E res -ascii  
else
  eval(['load ' SDIR 'BnMat.mat']);
end

mm = mk;
II = mm - mk(1) + 1;
mk = mk(II); 
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
   hf=figure(10*PlotBn + 0);
	res = abs(BnPEST)*1e+3;
   II = find(ss(:,1)<1.98);
   pcolor(mk,ss(II,:),res(II,:)), hold on, shading interp
   contour(mk,ss(II,:),res(II,:),'k-'), hold on,
   save B1spectrum mk ss res
   %axis([min(mk) max(mk) 0 1])
   axis([-15 15 0 1])
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('sqrt(\psi_p)','FontSize',16,'FontWeight','Bold')
   title('|\xi^1| [mm]','FontSize',16,'FontWeight','Bold')
   colorbar,  colormap(hot)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
   %map = colormap; N = size(map,1); a = linspace(0,1,N)';  b = a.^0.5*ones(1,3); colormap(map.*b)

   %plot rational surfaces
   if CheckQ > 0
      mq = [ceil(min(q)*abs(n)):max(mk)];
      qq = mq/abs(n);
      [sq,qn] = MacFindX(s,q,qq);
      mq = qn*n;
      plot(mq,sq,'b+','LineWidth',3,'MarkerSize',8), hold on,
      save Qsurf mq sq
      mm2 = mq; II2=mm2-mk(1)+1;
     %mm3 = [mq(end)+1:mq(end)+6]; II3=mm3-mk(1)+1; 
     %mm3 = mq(2); II3=mm3-mk(1)+1;
        mm3 = mk; II3=mm3-mk(1)+1;
      %mm3 = [1:4 6:17]; II3=mm3-mk(1)+1;
     %II3 = find(mk==11);  mm3=mk(II3);
     %II3 = II2;
   end

   if CheckQ > 0
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
   hf=figure(10*PlotBn + 1);
   s2 = ss.^2;
   s0 = 0.9;
   JJ = find(s2>s0);
   pcolor(mk,s2(JJ),1e+3*abs(BnPEST(JJ,:))), hold on, shading interp
   contour(mk,s2(JJ),1e+3*abs(BnPEST(JJ,:)),'k-'), hold on,
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

   hf=figure(10*PlotBn + 2);
   plot(ss.^2,real(BnPEST(:,II3))*1e+3,LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('Re(\xi^1_m) [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq].^2,[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 

   hf=figure(10*PlotBn + 3);
   plot(ss.^2,imag(BnPEST(:,II3))*1e+3,LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('Im(\xi^1_m) [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq].^2,[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 

   hf=figure(10*PlotBn + 4);
   for ks=1:length(II3)
   %hs=subplot(2,2,ks,'align');
   plot(ss.^2,abs(BnPEST(:,II3(ks)))*1e+3,LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('|\xi^1_m| [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   %set(hs,'FontSize',16,'FontWeight','Bold')

   IKINK = find(ss.^2<0.5);
   IPEEL = find(ss.^2>0.8);
   FIN_RES_Akink = max(max(abs(BnPEST(IKINK,:))))*1e+3;
   FIN_RES_Apeel = max(max(abs(BnPEST(IPEEL,:))))*1e+3;
   %IKINK = find(mk<5&mk>0);
   %IPEEL = find(mk>10);
   %FIN_RES_Akink = max(max(abs(BnPEST(:,IKINK))))*1e+3;
   %FIN_RES_Apeel = max(max(abs(BnPEST(:,IPEEL))))*1e+3;
   end

   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      for ks=1:length(II2)
      sqn = sq(mm2(ks)-mq(1)+1);
      plot([sqn.^2; sqn.^2],[a(3)*ones(size(sqn)); a(4)*ones(size(sqn))],'k--'), hold on,
      end   
   end 

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

      hf=figure(10*PlotBn + 11);
      plot(sq,abs(res(1,:)+res(2,:)*i)*1e+3,LSS,'LineWidth',2,'Color',LSC,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
      xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
      ylabel('|\xi^1_{res}| [mm]','FontSize',16,'FontWeight','Bold')
      ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
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

   if CheckQ > 0
   hf=figure(10*PlotBn + 5);
   JJ = mq - mk(1) + 1;
   plot(mk(JJ),1e+3*max(abs(BnPEST(:,JJ)),[],1),LSS,'Color',LSC,'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
   plot(mk,1e+3*max(abs(BnPEST),[],1),LSS,'Color',LSC,'LineWidth',2,'MarkerSize',7), hold on,
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('max|\xi^1_m(\psi_p)| [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   end

%  total Bn along plasma surface
   hf=figure(10*PlotBn + 7);
   [x,I]=sort(Tg); 
   plot(Tg(I)*180/pi,1e+3*real(BnEDGE(I)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('geometric \theta','FontSize',16,'FontWeight','Bold')
   ylabel('Re[\xi_n] [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

%  total Bn along plasma surface
   hf=figure(10*PlotBn + 8);
   [x,I]=sort(Tg); 
   plot(Tg(I)*180/pi,1e+3*abs(BnEDGE(I)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('geometric \theta','FontSize',16,'FontWeight','Bold')
   ylabel('|\xi_n| [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

   X = Tg(I)*180/pi;  Y = 1e+3*abs(BnEDGE(I));
   [Ymin,IKINK] = min(abs(X-0));
   IPEEL = find(X<-40);
   FIN_RES_DX = max(Y(IPEEL));
   FIN_RES_DM = Y(IKINK);

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

end

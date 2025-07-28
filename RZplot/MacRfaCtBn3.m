%function MacRfaCtBn2
%% The 2nd (fast) method to
%% transform from B1M in MARS coordinate to BnM=B1M in HAMADA coordinate
%% and plot the BnM poloidal harmonics
%% before running this subroutine, make sure that 
%% MacMain* is a procedure, not a function
%%   1) RMZM_F & BPLASMA are used in the MacMain* subroutine
%%   2) SDIR defined here agrees with those in the MacMain* subroutine 
%%   3) files RMZM_F_EQAC, BPLASMA(EQAC) exist in SDIR
%%   4) comment out Mac.RunB in MacMain*
%%   5) FEEDI is explicitely assumed to be =1 in MacMain*
%% NB: BnPEST is actually B1 in Hamada coordinates. 

global Mac 
global SDIR 

%SDIR = '/.automount/funsrv1/root/home/yliu/D3D_MATT/';        facn=pi/2;        SMAIN='D3D134234';
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST020333/New/';  facn=pi/2*1.7476; SMAIN='MAST020333';
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST_ALL/';  facn=pi/2*1.7476; SMAIN='MAST_ALL';
%SDIR = '/.automount/funsrv1/root/home/yliu/ANATOR/';  facn=1.0; SMAIN='ANATOR';
%SDIR = '/.automount/funsrv1/root/home/yliu/Rfx/ResEQ4b/';  facn=1.0; SMAIN='Rfx';
%SDIR = '/.automount/funsrv1/root/home/yliu/D3D117327/';       facn=pi/2;        SMAIN='D3D117327';
%SDIR = '/.automount/funsrv1/root/home/yliu/Scen2_V01/';  facn=1.0; SMAIN='ITER_15MA';
%SDIR = '/.automount/funsrv1/root/home/yliu/D3D139571/Rfa/';  facn=1.0; SMAIN='D3D139571';
SDIR = '/.automount/funsrv1/root/home/yliu/JET77329/';  facn=1.0; SMAIN='JET77329';
%SDIR = '/.automount/funsrv1/root/home/yliu/MAST22422/';  facn=1.0; SMAIN='MAST22422';


Mac.RunJ  = 0;
Mac.RunV  = 0;
n         = 1;
CheckPlot = 0;
CheckQ    = 1;
CheckVac  = 0;
CheckErgos= 0;
PlotBn    = 5;
PlotQ     = 3;
PlotPhi   = 0;
LSS       = '-o'; 
%LSC       = [(kclr-1)/(Nclr-1) (kclr-1)*(Nclr-kclr)*4/Nclr^2 (Nclr-kclr)/(Nclr-1)]
LSC       = [0 0 1];

mk = -29:29;
mi = mk - mk(1) + 1;
  
if CheckVac == 1
  eval(['load ' SDIR 'BnMat_VAC.mat']);
  BMnPEST_VAC  = BMnPEST;
  BnEDGE_VAC   = BnEDGE;
end

if CheckPlot == 0
  eval(['!cp ' SDIR 'RMZM_F_EQAC ' SDIR 'RMZM_F']);
  Mac.RunB = 1;
  eval(['MacMain' SMAIN]);
  II = 1:Mac.Ns1;
  BnEQAC = B1(II,:);
  BnEDGE = BnEQAC(end,:);
  ss     = Mac.s(II);

  BnEQAC(2:end,:) = BnEQAC(2:end,:)./jacobian(2:Mac.Ns1,:);
  
  MacInputPROFEQ([SDIR 'PROFEQ_PEST']);
  [t_hamada,f_hamada,tc_hamada] = MacGetHamada(R,jacobian);  
  BMnPEST = MacFunHamada(BnEQAC,t_hamada,f_hamada,tc_hamada,mk);

  save BnMat ss mk BMnPEST BnEDGE
  eval(['!mv ' 'BnMat.mat ' SDIR]);
else
  eval(['load ' SDIR 'BnMat.mat']);
end

mm = mk;
mm2= 2:4;
II = mm - mk(1) + 1;
mk = mk(II); 
BnPEST  = BMnPEST(:,II)/facn;
BnEDGE  = BnEDGE/facn;
BnPEST(1,:) = BnPEST(2,:);
if CheckVac == 1
   BnPEST_VAC  = BMnPEST_VAC(:,II)/facn;
   BnPEST_VAC(1,:) = BnPEST_VAC(2,:);
   BnEDGE_VAC  = BnEDGE_VAC/facn;
end
II2 = mm2 - mk(1) + 1;

%II = find(abs(BnPESTN)>0.3); BnPESTN(II)=0.3;

if CheckQ > 0
   s = ss(1:Mac.Ns1);
   q = Mac.q;
   dpsi = Mac.dpsi;

   % normalisation according to ERGOS
   if CheckErgos>0
      for k=1:size(BnPEST,2)
          BnPEST(2:end,k) = BnPEST(2:end,k)*2*pi./q(2:end)./dpsi(2:end);
      end
      if CheckVac>0
         for k=1:size(BnPEST,2)
             BnPEST_VAC(2:end,k) = BnPEST_VAC(2:end,k)*2*pi./q(2:end)./dpsi(2:end);
         end
      end
   end
end


if PlotBn > 0
   hf=figure(10*PlotBn + 0);
   pcolor(mk,ss,abs(BnPEST)), hold on, shading interp
   contour(mk,ss,abs(BnPEST),'k-'), hold on,
   axis([min(mk) max(mk) 0 1])
   xlabel('m','FontSize',14)
   ylabel('sqrt(\psi_p)','FontSize',14)
   title('|b^1_{tot}| [G/kAt] in PEST CS','FontSize',14)
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
		mm2 = mq; II2=mm2-mk(1)+1;
   end

   if CheckQ > 0
      hf=figure(10*PlotQ + 0);
      plot(s,q,LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
      xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
      ylabel('q','FontSize',16,'FontWeight','Bold')
      ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
      a = axis;
      for k=1:length(sq)
          plot([0 1],[mq(k) mq(k)],'k--'), hold on,
          plot([sq(k) sq(k)],[a(3) a(4)],'k--'), hold on,
      end
   end

   %plot Bn near the plasma boundary
	hf=figure(10*PlotBn + 1);
   s2 = ss.^2;
   s0 = 0.9;
   JJ = find(s2>s0);
   pcolor(mk,s2(JJ),abs(BnPEST(JJ,:))), hold on, shading interp
   contour(mk,s2(JJ),abs(BnPEST(JJ,:)),'k-'), hold on,
   axis([min(mk) max(mk) s0 1])
   xlabel('m','FontSize',14)
   ylabel('\psi_p','FontSize',14)
   title('|b^1_{tot}| [G/kAt] in PEST CS','FontSize',12)
   colorbar,  colormap(hot)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
   if CheckQ > 0
      k = find(sq>s0);
      plot(mq(k),sq(k).^2,'b+','LineWidth',3,'MarkerSize',8), hold on,
   end

   hf=figure(10*PlotBn + 2);
	plot(ss.^2,real(BnPEST(:,II2)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('Re(b^1_m)','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq].^2,[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 
   if CheckVac==1
   plot(ss.^2,real(BnPEST_VAC(:,II2)),'b--','LineWidth',2), hold on,
   end

   hf=figure(10*PlotBn + 3);
   plot(ss.^2,imag(BnPEST(:,II2)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('Im(b^1_m)','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckVac==1
   plot(ss.^2,imag(BnPEST_VAC(:,II2)),'b--','LineWidth',2), hold on,
   end
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq].^2,[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 

   hf=figure(10*PlotBn + 4);
   for ks=1:length(II2)
   %hs=subplot(2,2,ks,'align');
   plot(ss.^2,abs(BnPEST(:,II2(ks))),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('|b^1_m| (Gauss/kA)','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   %set(hs,'FontSize',16,'FontWeight','Bold')
   if CheckVac==1
    plot(ss.^2,abs(BnPEST_VAC(:,II2(ks))),'k--','LineWidth',2), hold on,
   end
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
      res = sq; resv = sq; peaka=sq;
      for k=1:length(sq)
          [X,ISQ] = min(abs(ss-sq(k)));
          res(k)  = abs(BnPEST(ISQ,mq(k)-mk(1)+1));
          peaka(k) = max(abs(BnPEST(:,mq(k)-mk(1)+1)));
          if CheckVac==1
             resv(k) = abs(BnPEST_VAC(ISQ,mq(k)-mk(1)+1));
          end
      end
      BnPEST_RS = [mq; sq; peaka; resv; res]'
      dpsiq     = spline(s,dpsi,sq);

      hf=figure(10*PlotBn + 11);
      plot(sq,res,LSS,'LineWidth',2,'Color',LSC,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
      if CheckVac==1
         plot(sq,resv,'--o','LineWidth',2,'Color','k','MarkerSize',7,'MarkerFaceColor','k'), hold on,
      end
      xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
      ylabel('|b^1_{res}| [G/kAt]','FontSize',16,'FontWeight','Bold')
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

   hf=figure(10*PlotBn + 5);
   JJ = mq - mk(1) + 1;
   plot(mk(JJ),max(abs(BnPEST(:,JJ)),[],1),LSS,'Color',LSC,'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
   plot(mk,max(abs(BnPEST),[],1),LSS,'Color',LSC,'LineWidth',2,'MarkerSize',7), hold on,
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('max|b^1_m(\psi_p)| [G/kAt]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckVac==1
   plot(mk(JJ),max(abs(BnPEST_VAC(:,JJ)),[],1),'--o','Color','k','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','k'), hold on,
   plot(mk,max(abs(BnPEST_VAC),[],1),'--o','Color','k','LineWidth',2,'MarkerSize',7), hold on,
   end

   %compare with vacuum field
   if CheckVac==1
	hf=figure(10*PlotBn + 6);
   plot(mk(JJ),max(abs(BnPEST(:,JJ)),[],1)./max(abs(BnPEST_VAC(:,JJ)),[],1),LSS,'Color',LSC,'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
   plot(mk,max(abs(BnPEST),[],1)./max(abs(BnPEST_VAC),[],1),LSS,'Color',LSC,'LineWidth',2,'MarkerSize',7), hold on,
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('max|b^{1tot}_m(\psi_p)|/max|b^{1vac}_m(\psi_p)|','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   end

   %pure plasma response
   if CheckVac==1
   hf=figure(10*PlotBn + 9);
   plot(ss.^2,abs(BnPEST(:,II2)-BnPEST_VAC(:,II2)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('|b^1_m(tot)-b^1_m(vac)|','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq.^2; sq.^2],[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 
   end

   %pure plasma response
   if CheckVac==1 & 1==0
   figure(10*PlotBn + 10)
   pcolor(mk,ss,abs(BnPEST-BnPEST_VAC)), hold on, shading interp
   contour(mk,ss,abs(BnPEST-BnPEST_VAC),'k-'), hold on,
   axis([min(mk) max(mk) 0 1])
   xlabel('m','FontSize',14)
   ylabel('sqrt(\psi_p)','FontSize',14)
   title('|b^1_{tot}-b^1_{vac}| [Gauss/kAt] in PEST CS','FontSize',12)
   colorbar,  colormap(hot)
   if CheckQ > 0
      plot(mq,sq,'b+','LineWidth',3,'MarkerSize',8), hold on,
   end
   end

   % field spectrum near the plasma surface
   IN = isnan(abs(BnPEST(:,1))); 
   II = find(IN==0);
   %II = find(ss<=0.99);
   IE = II(end);
   JJ = mq - mk(1) + 1;
   hf=figure(10*PlotBn + 13);
   plot(mk,abs(BnPEST(IE,:)),LSS,'LineWidth',1,'MarkerSize',7,'Color',LSC), hold on,
   plot(mk(JJ),abs(BnPEST(IE,JJ)),LSS,'Color',LSC,'LineWidth',1,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('|b^1_m(\psi_p=1)| [G/kAt]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   Ab1Rat = abs(BnPEST(IE,JJ(end)))/max(abs(BnPEST(IE,:)));
   if CheckVac==1
   plot(mk,abs(BnPEST_VAC(IE,:)),'b--o','LineWidth',2,'MarkerSize',8), hold on,
   plot(mk(JJ),abs(BnPEST_VAC(IE,JJ)),'bo','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','b'), hold on,
   %plot(mk,abs(BnPEST(IE,:)-BnPEST_VAC(IE,:)),'b--o','LineWidth',2,'MarkerSize',8), hold on,
   %plot(mk(JJ),abs(BnPEST(IE,JJ)-BnPEST_VAC(IE,JJ)),'bo','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','b'), hold on,
   Ab1RatV = abs(BnPEST_VAC(IE,JJ(end)))/max(abs(BnPEST_VAC(IE,:)));
   Ab1Rat = [Ab1RatV Ab1Rat];
   end
   res_b1res_ration = [Ab1Rat mk(JJ(end))]

   % displacement variation over phi angle
   if PlotPhi > 0
   hf=figure(10*PlotPhi + 0);
   x = Mac.phi*180/pi;
   y = real(BnEDGE_VAC(IXL)*exp(-i*n*Mac.phi)); y=y/max(y);
   plot(x,y,'--','Color','k','LineWidth',1), hold on,
   end

end
  

  



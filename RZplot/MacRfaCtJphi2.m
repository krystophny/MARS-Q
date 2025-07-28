%function MacRfaCtBn2
%% The 2nd (fast) method to
%% transform from B1M in MARS coordinate to BnM in PEST coordinate
%% and plot the BnM poloidal harmonics
%% before running this subroutine, make sure that 
%% MacMain* is a procedure, not a function
%%   1) RMZM_F & BPLASMA are used in the MacMain* subroutine
%%   2) SDIR defined here agrees with those in the MacMain* subroutine 
%%   3) files RMZM_F_EQAC, RMZM_F_PEST, BPLASMA(EQAC) exist in SDIR
%%   4) comment out Mac.RunB in MacMain*
%%   5) FEEDI is explicitely assumed to be =1 in MacMain*
%% NB: BnPEST is actually B1PEST=J*(b\cdot\nabla s). Only B1PEST has the propoerty 
%% that of vanishing amplitude at the rational surfaces in ideal MHD

global Mac 
global SDIR FEEDI

%SDIR = '/.automount/funsrv1/root/home/yliu/D3D_MATT/';        facn=pi/2;        SMAIN='D3D134234';
SDIR = '/.automount/funsrv1/root/home/yliu/MAST020333/New/';  facn=pi/2*1.7476; SMAIN='MAST020333';
%SDIR = '/.automount/funsrv1/root/home/yliu/ANATOR/';  facn=1.0; SMAIN='ANATOR';
%SDIR = '/.automount/funsrv1/root/home/yliu/Rfx/ResEQ4b/';  facn=1.0; SMAIN='Rfx';
%SDIR = '/.automount/funsrv1/root/home/yliu/D3D117327/';       facn=pi/2;        SMAIN='D3D117327';
%SDIR = '/.automount/funsrv1/root/home/yliu/Scen2_V01/';  facn=1.0; SMAIN='ITER_15MA';
Mac.RunB  = 0;
n         = 3;
CheckPlot = 1;
CheckQ    = 1;

PlotBn    = 7;

mk = -29:29;
mi = mk - mk(1) + 1;
  
if CheckPlot == 0
  eval(['!cp ' SDIR 'RMZM_F_EQAC ' SDIR 'RMZM_F']);
  Mac.RunJ = 1;
  eval(['MacMain' SMAIN]);
  II = find(Mac.s<=1.0&Mac.s>=0.0);
  BnEQAC = Jphi(II,:);
  R_EQAC = R(II,:);
  Z_EQAC = Z(II,:);

  %compute Bn vs. geometrical angle theta at the plasma surface
  Rs = R_EQAC(end,:);
  Zs = Z_EQAC(end,:);
  Rc = (min(Rs)+max(Rs))/2;
  Zc = (min(Zs)+max(Zs))/2;
  Tg = atan2(Zs-Zc,Rs-Rc);
  BnEDGE = BnEQAC(end,:);

  eval(['!cp ' SDIR 'RMZM_F_PEST ' SDIR 'RMZM_F']);
  Mac.RunJ = 0;
  eval(['MacMain' SMAIN]);
  ss     = Mac.s(II);
  R_PEST = R(II,:);
  Z_PEST = Z(II,:);
  %G22_PEST  = dRdchi(II,:).^2 + dZdchi(II,:).^2;  G22_PEST(1,:) = G22_PEST(2,:);

  BnPEST  = griddata(R_EQAC,Z_EQAC,BnEQAC,R_PEST,Z_PEST);
  %BnPEST = BnPEST.*sqrt(G22_PEST).*R_PEST;  %B1PEST
  expmchi = exp(-Mac.chi'*mk*i);
  BMnPEST = BnPEST*expmchi*(Mac.chi(2)-Mac.chi(1))/2/pi;

  save JphiMat ss mk Tg BMnPEST BnEDGE
  eval(['!mv ' 'JphiMat.mat ' SDIR]);
else
  eval(['load ' SDIR 'JphiMat.mat']);
end

mm = -29:29;
mm2= [5];
II = mm - mk(1) + 1;
mk = mk(II); 
BnPEST  = BMnPEST(:,II)/facn;
BnEDGE  = BnEDGE/facn;
BnPEST(1,:) = BnPEST(2,:);
II2 = mm2 - mk(1) + 1;

%II = find(abs(BnPESTN)>0.3); BnPESTN(II)=0.3;

if PlotBn > 0
   hf=figure(10*PlotBn + 0);
   pcolor(mk,ss,abs(BnPEST)), hold on, shading interp
   contour(mk,ss,abs(BnPEST),'k-'), hold on,
   axis([min(mk) max(mk) 0 1])
   xlabel('m','FontSize',14)
   ylabel('sqrt(\psi_p)','FontSize',14)
   title('|j_\phi| [kA/kA] in PEST CS','FontSize',14)
   colorbar,  colormap(hot)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
   %map = colormap; N = size(map,1); a = linspace(0,1,N)';  b = a.^0.5*ones(1,3); colormap(map.*b)

   %plot rational surfaces
   if CheckQ > 0
      dataq = load([SDIR 'PROFEQ_PEST']);
      %dataq = load([SDIR 'PROFQ_PEST']);
      s = dataq(:,1);
      q = dataq(:,2);
      mq = [ceil(min(q)*abs(n)):max(mk)];
      qq = mq/abs(n);
      [sq,qn] = MacFindX(s,q,qq);
      mq = qn*abs(n);
      plot(mq,sq,'b+','LineWidth',3,'MarkerSize',8), hold on,

      IQ = [4 14];
      mq = mq(IQ); sq = sq(IQ);
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
   title('|j_\phi| [kA/kA] in PEST CS','FontSize',12)
   colorbar,  colormap(hot)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
   if CheckQ > 0
      k = find(sq>s0);
      plot(mq(k),sq(k).^2,'b+','LineWidth',3,'MarkerSize',8), hold on,
   end

   hf=figure(10*PlotBn + 2);
	plot(ss.^2,real(BnPEST(:,II2)),'r-','LineWidth',2), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('Re(j_{\phi,m}) [kA/kA]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq].^2,[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 

   hf=figure(10*PlotBn + 3);
   plot(ss.^2,imag(BnPEST(:,II2)),'r-','LineWidth',2), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('Im(j_{\phi,m}) [kA/kA]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq].^2,[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 

   hf=figure(10*PlotBn + 4);
   plot(ss.^2,abs(BnPEST(:,II2)),'r-','LineWidth',2), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('|j_{\phi,m}| (kA/kA)','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq.^2; sq.^2],[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 
	%find abs(BnPEST) at rational surfaces
   if CheckQ > 0
      res = sq;
      for k=1:length(sq)
          [X,ISQ] = min(abs(ss-sq(k)));
          res(k)  = abs(BnPEST(ISQ,mq(k)-mk(1)+1));
      end
      BnPEST_RS = [mq; sq.^2; res]
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
   plot(mk,max(abs(BnPEST),[],1),'r-o','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','r'), hold on,
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('max|j_{\phi,m}(\psi_p)| [kA/kA]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

%  Bn along plasma surface
   hf=figure(10*PlotBn + 7);
   plot(Tg(2:end),real(BnEDGE(2:end)),'r-','LineWidth',2), hold on,
   xlabel('\theta','FontSize',16,'FontWeight','Bold')
   ylabel('Re[j_\phi] [kA/kA]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

end
  

  



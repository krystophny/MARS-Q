function [BresQ3]=MacRfaCtRMPthreshold;  
%% The 2nd (fast) method to
%% transform from B1M in MARS coordinate to BnM in PEST coordinate
%% and plot the BnM poloidal harmonics
%% before running this subroutine, make sure that 
%% MacMain* is a procedure, not a function
%%   1) RMZM_F & BPLASMA are used in the MacMain* subroutine
%%   2) SDIR defined here agrees with those in the MacMain* subroutine 
%%   3) files RMZM_F_EQAC, RMZM_F_PEST, BPLASMA_EQAC exist in SDIR
%%   4) comment out Mac.RunB in MacMain*
%%   5) FEEDI is explicitely assumed to be =1 in MacMain*
%% NB: BnPEST is actually B1PEST=J*(b\cdot\nabla s). Only B1PEST has the property 
%%     of vanishing amplitude at the rational surfaces in ideal MHD

global Mac 
global SDIR SDIR5 k_plot
global SDIRW

SMAIN='ALL_B';

facn      = 1.0;
Mac.RunJ  = 0;
Mac.RunV  = 0;
CheckQ    = 1;
CheckErgos= 0;  %=1: ERGOS definition; =2: SURFMN definition; =3: remove q-factor from b1
PlotBn    = 0;
PlotB1E   = 0;
PlotQ     = 0;
PlotPhi   = 0;
LSS       = '-o'; 
LSC       = [0 1 0];

if k_plot==1, PlotBn=2; PlotQ=3; end

if 1==1 
mk = -40:40;
mi = mk - mk(1) + 1;
  
copyfile([SDIR 'RMZM_F.OUT'],[SDIR 'RMZM_F_EQAC.OUT'],'f');
Mac.RunB = 1;
eval(['MacMain' SMAIN]);
II = 1:Mac.Ns1;
BnEQAC = Bn(II,:);
R_EQAC = R(II,:);
Z_EQAC = Z(II,:);

n = abs(Mac.n);

%compute Bn vs. geometrical angle theta at the plasma surface
Rs = R_EQAC(Mac.Ns1,:);
Zs = Z_EQAC(Mac.Ns1,:);
Rc = (min(Rs)+max(Rs))/2;
Zc = (min(Zs)+max(Zs))/2;
Tg = atan2(Zs-Zc,Rs-Rc);
[Y,IXU] = max(Zs); 
[Y,IXL] = min(Zs); 
[Y,IXM] = max(Rs); 
BnEDGE = BnEQAC(Mac.Ns1,:);

copyfile([SDIR 'RMZM_F_PEST.OUT'],[SDIR 'RMZM_F.OUT'],'f');
Mac.RunB = 0;
eval(['MacMain' SMAIN]);
ss     = Mac.s(II);
R_PEST = R(II,:);
Z_PEST = Z(II,:);
G22_PEST  = dRdchi(II,:).^2 + dZdchi(II,:).^2;  G22_PEST(1,:) = G22_PEST(2,:);

BnPEST  = griddata(R_EQAC,Z_EQAC,BnEQAC,R_PEST,Z_PEST);
BnPEST = BnPEST.*sqrt(G22_PEST).*R_PEST;  %B1PEST

copyfile([SDIR 'RMZM_F_EQAC.OUT'],[SDIR 'RMZM_F.OUT'],'f');

ss = Mac.s(1:Mac.Ns1);
BnPEST = BnPEST(1:Mac.Ns1,:);
BnPEST2 = BnPEST./jacobian(1:Mac.Ns1,:);
BnPEST2(1,:) = BnPEST2(2,:); 

expmchi = exp(-Mac.chi'*mk*i);
BMnPEST = BnPEST*expmchi*(Mac.chi(2)-Mac.chi(1))/2/pi;

B1E = BnPEST(end,:); B1E=B1E(:);

end


%patch possible NaN's of BMnPEST near the plasma boundary
if 1==1
for k=1:size(BMnPEST,2)
    INOR = find(isnan(BMnPEST(:,k))==0);
    if length(INOR) < size(BMnPEST,1)
       II = INOR(end); 
       BMnPEST(II+1:end,k) = BMnPEST(II,k);
    end
end
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
   if CheckErgos==1,  fac_ERGOS = q.*dpsi/2; fac_ERGOS(1)=fac_ERGOS(2); end
   if CheckErgos==2,  fac_ERGOS = q./dpsi; end
   if CheckErgos==3,  fac_ERGOS = q; end
   if CheckErgos>0
      for k=1:size(BnPEST,2)
          BnPEST(1:end,k) = BnPEST(1:end,k)./fac_ERGOS(1:end);
      end
   end
end


%plot rational surfaces
if CheckQ > 0
   mq = [ceil(min(q)*abs(n)):max(mk)];
   qq = mq/abs(n);
   [sq,qn] = MacFindX(s,q,qq);
   mq = qn*n;
   mm2 = mq; II2=mm2-mk(1)+1;
   II3 = II2;
end

%find abs(BnPEST) at rational surfaces
if CheckQ > 0
res = [sq; sq]; resv = sq; peaka=sq;
for k=1:length(sq)
    [X,ISQ] = min(abs(ss-sq(k)));
    res(:,k)  = [real(BnPEST(ISQ,mq(k)-mk(1)+1)); imag(BnPEST(ISQ,mq(k)-mk(1)+1))];
    peaka(k)= max(abs(BnPEST(:,mq(k)-mk(1)+1)));
end
format short e
BnPEST_RS  = transpose([mq; sq; peaka; res; abs(res(1,:)+res(2,:)*i)])
dpsiq      = spline(s,dpsi,sq);

CRIT_b1res = [res(1,end) res(2,end)];
save([SDIR5 'CRIT_b1res'],'CRIT_b1res','-ascii');

a = [0 0]; NA=0;
for k=1:length(sq)
    if sq(k)>0.975 | k==length(sq), a=a+[res(1,k) res(2,k)]; NA=NA+1; end
end
CRIT_b1res_mean = a/NA;
save([SDIR5 'CRIT_b1res_mean'],'CRIT_b1res_mean','-ascii');

if PlotBn > 0
   if 1==0
   hf=figure(10*PlotBn + 0);
   resBn = abs(BnPEST);
   II = find(ss(:,1)<1.00);
   pcolor(mk,ss(II,:),resBn(II,:)),  colorbar, hold on, shading interp
  %caxis([-0.1 12.1])
 %contour(mk,ss(II,:),resBn(II,:),'k-'), hold on,
   save B1spectrum mk ss resBn
   %axis([min(mk) max(mk) 0 1])
   axis([-20 20 0 1]), 
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('sqrt(\psi_p)','FontSize',14)
   title('|b^1|x10^4','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',12)
   colormap(hot)
   %map = colormap; N = size(map,1); a = linspace(0,1,N)';  b = a.^0.5*ones(1,3); colormap(map.*b)
   end

   if CheckQ > 0
      hf=figure(10*PlotQ + 0);
      plot(s.^2,q,LSS(1:end-1),'LineWidth',3,'Color',LSC), hold on,
      xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
      ylabel('q','FontSize',18,'FontWeight','Bold')
      ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
      a = axis;
      for k=1:length(sq)
          plot([0 1],[mq(k) mq(k)]/n,'k--'), hold on,
          plot([sq(k)^2 sq(k)^2],[a(3) a(4)],'k--'), hold on,
      end
      saveas(hf,'profq.fig')

      if 1==0
      dw=load([SDIR 'PROFROT.IN']);
%      de=load([SDIR 'PROFWE.IN']);
      hf=figure(10*PlotQ + 1);
      plot(dw(2:end,1).^2,dw(2:end,2),'r-','LineWidth',3), hold on,
%      plot(de(2:end,1).^2,de(2:end,2),'b--','LineWidth',3), hold on,
      xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
      ylabel('rotation frequency [rad/s]','FontSize',18,'FontWeight','Bold')
      ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
      legend('\Omega','\Omega_{ExB}')
      a = axis;
      for k=1:length(sq)
          plot([0 1],[mq(k) mq(k)]/n,'k--'), hold on,
          plot([sq(k)^2 sq(k)^2],[a(3) a(4)],'k--'), hold on,
      end
      saveas(hf,'profw.fig')
      end
  end
   
   %plot Bn near the plasma boundary
   hf=figure(10*PlotBn + 1);
   s2 = ss.^1;
   s0 = 0.0;
   JJ = find(s2>s0);
   pcolor(mk,s2(JJ),abs(BnPEST(JJ,:))), hold on, shading interp
   %axis([-10 10 s0 1])
   %save IterCase1Btfte2.mat mk ss BnPEST
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   title('|b^1|x10^4','FontSize',16,'FontWeight','Bold')
   %title('|b_n| [Gauss]','FontSize',16,'FontWeight','Bold')
%   colorbar,  colormap(load('colorbar.txt')')
   colorbar,  colormap(hot)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
   if CheckQ > 0
      k = find(sq>s0);
      plot(mq(k),sq(k).^1,'b+','LineWidth',3,'MarkerSize',9), hold on,
   end
   saveas(hf,'b1m2d.fig')

   if 1==0
   hf=figure(10*PlotBn + 2);
   plot(ss.^2,real(BnPEST(:,II3)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('Re(b^1_m)','FontSize',16,'FontWeight','Bold')
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
   ylabel('Im(b^1_m)','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq].^2,[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 
   end

   hf=figure(10*PlotBn + 4);
   for ks=1:length(II3)
   %hs=subplot(2,2,ks,'align');
   plot(ss.^2,abs(BnPEST(:,II3(ks))),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('|b^1_m| x 10^4','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   %set(hs,'FontSize',16,'FontWeight','Bold')
   end

   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      for ks=1:length(II2)
      %sqn = sq(round(mm2(ks)-mq(1)+1));
      sqn = sq(ks);
      plot([sqn.^2; sqn.^2],[a(3)*ones(size(sqn)); a(4)*ones(size(sqn))],'k--'), hold on,
      end   
   end 
   saveas(hf,'b1m1d.fig')

   if CheckQ > 0 
      hf=figure(10*PlotBn + 11);
      plot(sq.^2,abs(res(1,:)+res(2,:)*i),LSS,'LineWidth',2,'Color',LSC,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
      xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
      ylabel('|b^1_{res}|x10^4 ','FontSize',16,'FontWeight','Bold')
      ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

      FIN_RES_b1res = abs(res(1,end)+res(2,end)*i)*1e+4;
      %FIN_RES_b1res = abs(res(1,end-1)+res(2,end-1)*i)*1e+4;
      %FIN_RES_b1res = abs(res(1,end-2)+res(2,end-2)*i)*1e+4;
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

   if 1==0
   if CheckQ > 0
   hf=figure(10*PlotBn + 5);
   JJ = mq - mk(1) + 1;
   plot(mk(JJ),max(abs(BnPEST(:,JJ)),[],1),LSS,'Color',LSC,'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
   plot(mk,max(abs(BnPEST),[],1),LSS,'Color',LSC,'LineWidth',2,'MarkerSize',7), hold on,
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('max|b^1_m(\psi_p)| [G/kAt]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   end
   end

%  total Bn along plasma surface
   if 1==0
   hf=figure(10*PlotBn + 7);
   [x,I]=sort(Tg); 
   plot(Tg(I),real(BnEDGE(I)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('geometric \theta','FontSize',16,'FontWeight','Bold')
   ylabel('Re[b_n] [Gauss/kAt]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   end

%  total Bn along plasma surface
   if 1==0
   hf=figure(10*PlotBn + 8);
   [x,I]=sort(Tg); 
   plot(Tg(I),abs(BnEDGE(I)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('geometric \theta','FontSize',16,'FontWeight','Bold')
   ylabel('|b_n| [Gauss/kAt]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   end

   % field spectrum near the plasma surface
   if CheckQ > 0 & PlotB1E>0
   IN = isnan(abs(BnPEST(:,1))); 
   II = find(IN==0);
   %II = find(ss<=0.995);
   IE = II(end);
   JJ = mq - mk(1) + 1;
   hf=figure(10*PlotB1E + 1);
   plot(mk,real(BnPEST(IE,:)),LSS,'LineWidth',1,'MarkerSize',7,'Color',LSC), hold on,
   plot(mk(JJ),real(BnPEST(IE,JJ)),LSS,'Color',LSC,'LineWidth',1,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('|b^1_m(\psi_p=1)| [G/kAt]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   end

   %  plot B1 at plasma surface, along PEST poloidal angle 
   if PlotB1E>0
   hf=figure(10*PlotB1E + 2);
   plot(Mac.chi,real(B1E),'-',Mac.chi,imag(B1E),'--','LineWidth',2,'Color',LSC), hold on,
   xlabel('\theta_{PEST}','FontSize',16,'FontWeight','Bold')
   ylabel('b^1 [Gauss/kAt]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   legend('real','imaginary')
   end

   %reconstruct (approximate) plasma displacement from Q1
    if CheckQ > 0 & PlotB1E>0
    expmchi = exp(mk'*Mac.chi*i);
    X1E = (BnPEST(IE,:)./(mk-n*q(IE)))*expmchi;
   hf=figure(10*PlotB1E + 3);
   [x,I]=sort(Tg); 
   plot(Tg(I),abs(X1E(I)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('geometric \theta','FontSize',16,'FontWeight','Bold')
   ylabel('\xi^1 [a.u.]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   end

%plot b1res in (R_PEST,Z_PEST) plane
   if 1==0
   hf=figure(10*PlotBn + 12);
	res = abs(BnPEST2)*1e+4;
   II = 1:Mac.Ns1;
   pcolor(R_PEST(II,:)*Mac.R0EXP,Z_PEST(II,:)*Mac.R0EXP,res(II,:)),  colorbar, hold on, shading interp
   xlabel('R [m]','FontSize',16,'FontWeight','Bold')
   ylabel('Z [m]','FontSize',16,'FontWeight','Bold')
   title('|b^1|x10^4','FontSize',14,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   %colormap(hot)
   axis equal
   end

end
end

function MacACADpest(SDIR2,SCOIL,SFLD,SEXB)
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
%% NB: BnPEST is actually B1PEST=J*(b\cdot\nabla s). Only B1PEST has the property 
%%     of vanishing amplitude at the rational surfaces in ideal MHD

global Mac Acad
global SDIR 

SDIR = SDIR2;

facn      = 1.0;
Mac.RunJ  = 0;
Mac.RunV  = 0;
n         = Acad.n;
CheckQ    = 1;
CheckErgos= 0;  %=1: ERGOS definition; =2: SURFMN definition; =3: remove q-factor from b1
PlotBn    = 8;
PlotB1E   = 0;
PlotQ     = 9;
PlotPhi   = 0;
LSS       = '-o'; 
%LSC       = [(kclr-1)/(Nclr-1) (kclr-1)*(Nclr-kclr)*4/Nclr^2 (Nclr-kclr)/(Nclr-1)]
LSC       = [0 0 1];

mk = -20:20;
mi = mk - mk(1) + 1;
  
copyfile([SDIR 'RMZM_F_EQAC'],[SDIR 'RMZM_F.OUT'],'f');
Mac.RunB = 1;
MacMainD3D_ALL_B;
II = 1:Mac.Ns;
BnEQAC = Bn(II,:);
R_EQAC = R(II,:);
Z_EQAC = Z(II,:);

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

copyfile([SDIR 'RMZM_F_PEST'],[SDIR 'RMZM_F.OUT'],'f');
Mac.RunB = 0;
MacMainD3D_ALL_B;
ss     = Mac.s(II);
R_PEST = R(II,:);
Z_PEST = Z(II,:);
G22_PEST  = dRdchi(II,:).^2 + dZdchi(II,:).^2;  G22_PEST(1,:) = G22_PEST(2,:);

BnPEST  = griddata(R_EQAC,Z_EQAC,BnEQAC,R_PEST,Z_PEST);
BnPEST = BnPEST.*sqrt(G22_PEST).*R_PEST;  %B1PEST

ss = Mac.s(1:Mac.Ns1);
BnPEST = BnPEST(1:Mac.Ns1,:);
BnPEST2 = BnPEST./jacobian(1:Mac.Ns1,:);
BnPEST2(1,:) = BnPEST2(2,:); 

expmchi = exp(-Mac.chi(:)*mk*i);
BMnPEST = BnPEST*expmchi*(Mac.chi(2)-Mac.chi(1))/2/pi;

if 1==1
  save BnMat ss mk Tg BMnPEST BnEDGE
  eval(['movefile ' 'BnMat.mat ' SDIR]);

  B1E = BnPEST(end,:); B1E=B1E(:);
  res = [Mac.chi(:) real(B1E) imag(B1E)];
  save FootPrintData_B1E res -ascii  
else
  eval(['load ' SDIR 'BnMat.mat']);
end

%patch possible NaN of BMnPEST near the plasma boundary
if 1==1
for k=1:size(BMnPEST,2)
    INOR = find(isnan(BMnPEST(:,k))==0);
    if length(INOR) < size(BMnPEST,1)
       %BMnPEST(:,k) = pchip(ss(INOR),BMnPEST(INOR,k),ss);
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

   Acad.Qprof =[s(:).^2 q(:)];

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


if PlotBn > 0
   %get rational surfaces
   if CheckQ > 0
      mq = [ceil(min(q)*abs(n)):max(mk)];
      qq = mq/abs(n);
      [sq,qn] = MacFindX(s,q,qq);
      %sq=[0.74958 0.94537 0.98740 0.99782].^0.5;
      mq = qn*n;
      %plot(mq,sq,'b+','LineWidth',3,'MarkerSize',8), hold on,
      save Qsurf mq sq
      mm2 = mq; II2=mm2-mk(1)+1;
     %mm3 = [mq(end)+1:mq(end)+6]; II3=mm3-mk(1)+1; 
     %mm3 = [2 3]; II3=mm3-mk(1)+1;
     II3 = II2;
   end

   if CheckQ > 0
      if 1==0
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
      end
  end
   
   %plot Bn near the plasma boundary
   hf=figure(10*PlotBn + 1);
   s2 = ss.^2;
   s0 = 0.5;
   JJ = find(s2>s0);
   pcolor(mk,s2(JJ),1e+4*abs(BnPEST(JJ,:))), hold on, shading interp
   %axis([-10 10 s0 1])
   %save IterCase1Btfte2.mat mk ss BnPEST
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ylabel('\psi_p','FontSize',16,'FontWeight','Bold')
   title('|b^1|x10^4','FontSize',16,'FontWeight','Bold')
   colorbar,  colormap(hot)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
   if CheckQ > 0
      k = find(sq>s0);
      plot(mq(k),sq(k).^1,'b+','LineWidth',3,'MarkerSize',8), hold on,
   end

   if 1==0
   hf=figure(10*PlotBn + 4);
   for ks=1:length(II3)
   %hs=subplot(2,2,ks,'align');
   plot(ss.^2,abs(BnPEST(:,II3(ks)))*1e+4,LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('|b^1_m| x 10^4','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   %set(hs,'FontSize',16,'FontWeight','Bold')
   end

   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      for ks=1:length(II2)
      sqn = sq(round(mm2(ks)-mq(1)+1));
      plot([sqn.^2; sqn.^2],[a(3)*ones(size(sqn)); a(4)*ones(size(sqn))],'k--'), hold on,
      end   
   end 
   saveas(hf,'b1m1d.fig')
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
      BnPEST_RS = [mq; sq; peaka; res; abs(res(1,:)+res(2,:)*i)]'
      dpsiq     = spline(s,dpsi,sq);

      hf=figure(10*PlotBn + 2);
      plot(sq.^2,abs(res(1,:)+res(2,:)*i)*1e+4,LSS,'LineWidth',2,'Color',LSC,'MarkerSize',7,'MarkerFaceColor',LSC), hold on,
      xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
      ylabel('|b^1_{res}|x10^4 ','FontSize',16,'FontWeight','Bold')
      ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

      FIN_RES_b1res = abs(res(1,end)+res(2,end)*i)*1e+4;
      %FIN_RES_b1res = abs(res(1,end-1)+res(2,end-1)*i)*1e+4;
      %FIN_RES_b1res = abs(res(1,end-2)+res(2,end-2)*i)*1e+4;

      x = sq.^2;
      y = abs(res(1,:)+res(2,:)*i)*1e+4;
      Acad.Bpest = [x(:) y(:)];
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

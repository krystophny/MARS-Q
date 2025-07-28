%function MacRfaCtBn1
%% The 1st (slow) method to 
%% transform from B1M in MARS coordinate to BnM in PEST coordinate
%% and plot the BnM poloidal harmonics
%% before running this subroutine, make sure that 
%% MacMain* is a procedure, not a function
%% 1) RMZM_F & BPLASMA are used in the MacMain* subroutine
%% 2) SDIR defined here agrees with those in the MacMain* subroutine 
%% 3) files RMZM_F_EQAC, RMZM_F_PEST, BPLASMA_EQAC exist in SDIR
%% 4) if BPLASMA_PEST also exist, set CheckPest = 1, otherwise = 0
%% 5) comment out FEEDI in MacMain*

global Mac 
global SDIR FEEDI

%SDIR = '/.automount/funsrv1/root/home/yliu/D3D_MATT/';   facn=pi/2;        FDI=[0.4735 0.3734]; SMAIN='D3D134234';
SDIR = '/.automount/funsrv1/root/home/yliu/MAST020333/New/';  facn=pi/2*1.7476; FDI=[1.0 1.0]; SMAIN='MAST020333';
%SDIR = '/.automount/funsrv1/root/home/yliu/D3D117327/';  facn=pi/2;        FDI=[0.4854 0.3645]; SMAIN='D3D117327';
CheckPest = 0;
CheckPlot = 1;
CheckQ    = 1;

mk = -29:29;
mi = mk - mk(1) + 1;
Mac.plot_Bn = 0;
Mac.plot_Bn2 = 5;
Mac.plot_Bn3 = 6;

if CheckPlot == 0
% define an array of radial points, where we do the coordinate transform
%rct = linspace(0,1,81); rct=rct(2:end); rct2 = rct;
%rct = 1.0;
rct  = Mac.s(2:Mac.Ns1)'; rct2=rct;

BnPEST  = zeros(length(mk),length(rct));
BnPESTN = BnPEST;
B1PEST  = BnPEST;
B1PESTN = BnPEST;
BnGEOM  = BnPEST;
BnGEOMN = BnPEST;

for k=1:length(rct)
  Mac.rbn = rct(k);
  eval(['!cp ' SDIR 'RMZM_F_EQAC ' SDIR 'RMZM_F']);
  eval(['!cp ' SDIR 'BPLASMA_EQAC ' SDIR 'BPLASMA']);
  FEEDI = FDI(1);
  Mac.SS = 'r-';
  eval(['MacMain' SMAIN]);
  [A1_EQAC,A2_EQAC,A3_EQAC,II] = MacGetCtMat(R,Z,dRdchi,dZdchi);
  A1_EQAC = A1_EQAC(mi,mi);
  A2_EQAC = A2_EQAC(mi,mi);
  A3_EQAC = A3_EQAC(mi,mi);
  B1M_EQAC = BM1(II,mi);

  eval(['!cp ' SDIR 'RMZM_F_PEST ' SDIR 'RMZM_F']);
  if CheckPest == 1
    eval(['!cp ' SDIR 'BPLASMA_PEST ' SDIR 'BPLASMA']);
  end
  FEEDI = FDI(2);
  Mac.SS = 'b--';
  eval(['MacMain' SMAIN]);
  [A1_PEST,A2_PEST,A3_PEST,II] = MacGetCtMat(R,Z,dRdchi,dZdchi);
  A1_PEST = A1_PEST(mi,mi);
  A2_PEST = A2_PEST(mi,mi);
  A3_PEST = A3_PEST(mi,mi);

  if CheckPest == 1
    B1M_PEST = BM1(II,mi);
    B1PEST(:,k) = B1M_PEST(:);
    BnPEST(:,k) = A2_PEST*B1M_PEST(:);
    BnGEOM(:,k) = A1_PEST*B1M_PEST(:);
  end

  BnGEOMN(:,k) = A1_EQAC*B1M_EQAC(:);
  B1PESTN(:,k) = A3_PEST*BnGEOMN(:,k);
  BnPESTN(:,k) = A2_PEST*B1PESTN(:,k);
 
  rct2(k) = Mac.s(II);
  res_status = [k rct2(k)]

  if Mac.plot_Bn > 0
    figure(10*Mac.plot_Bn + 0)
    plot(Mac.Mm(mi),abs(BnPESTN(:,k)),'r-o','LineWidth',1,'MarkerSize',7), hold on,
    xlabel('m','FontSize',16,'FontWeight','Bold')
    ylabel('|b_n(PEST)| [Gauss/kA]','FontSize',16,'FontWeight','Bold')
  
    if CheckPest == 1
       plot(Mac.Mm(mi),abs(BnPEST(:,k)),'b--s','LineWidth',1,'MarkerSize',7), hold on,
    end

    figure(10*Mac.plot_Bn + 1)
    plot(Mac.Mm(mi),abs(B1PESTN(:,k)),'r-o','LineWidth',1,'MarkerSize',7), hold on,
    xlabel('m','FontSize',16,'FontWeight','Bold')
    ylabel('|b_s(PEST)| [Gauss/kA]','FontSize',16,'FontWeight','Bold')
  
    if CheckPest == 1
       plot(Mac.Mm(mi),abs(B1PEST(:,k)),'b--s','LineWidth',1,'MarkerSize',7), hold on,
    end

    figure(10*Mac.plot_Bn + 2)
    plot(Mac.Mm(mi),abs(BnGEOMN(:,k)),'r-o','LineWidth',1,'MarkerSize',7), hold on,
    xlabel('m','FontSize',16,'FontWeight','Bold')
    ylabel('|b_n(GEOM)| [Gauss/kA]','FontSize',16,'FontWeight','Bold')
  
    if CheckPest == 1
       plot(Mac.Mm(mi),abs(BnGEOM(:,k)),'b--s','LineWidth',1,'MarkerSize',7), hold on,
    end
  end
end

rct = rct2;
end

if CheckPlot == 1
  eval(['load ' SDIR 'BnMat.mat']);
else
  save BnMat rct mk BnPEST BnPESTN B1PEST B1PESTN BnGEOM BnGEOMN
end

N0 = 0;
II = -mk(1)-14-N0:-mk(1)+16+N0;
mk = mk(II); 
BnPEST  = BnPEST(II,:);
BnPESTN = BnPESTN(II,:);
B1PEST  = B1PEST(II,:);
B1PESTN = B1PESTN(II,:);
BnGEOM  = BnGEOM(II,:);
BnGEOMN = BnGEOMN(II,:);

BnPEST  = BnPEST/facn;
BnPESTN = BnPESTN/facn;

%II = find(abs(BnPESTN)>0.3); BnPESTN(II)=0.3;

if Mac.plot_Bn2 > 0
   figure(10*Mac.plot_Bn2 + 0)
   pcolor(mk,rct,abs(BnPESTN)'), hold on, shading interp
	contour(mk,rct,abs(BnPESTN)','k-'), hold on,
   axis([min(mk) max(mk) 0 1])
   xlabel('m','FontSize',14)
   ylabel('sqrt(\psi_p)','FontSize',14)
   title('|b_n| [Gauss/kA] in PEST CS','FontSize',12)
   colorbar, map = colormap; N = size(map,1); 
   a = linspace(0,1,N)';  b = a.^0.5*ones(1,3); 
   %colormap(map.*b)
   colormap(hot)

   %plot q-profile
   if CheckQ > 0
      dataq = load([SDIR 'PROFEQ']);
      s = dataq(:,1);
      q = dataq(:,2);
      mq = [ceil(min(q)*abs(Mac.n)):max(mk)];
      qq = mq/abs(Mac.n);
		[sq,qn] = MacFindX(s,q,qq);
      mq = qn*abs(Mac.n);
      plot(mq,sq,'b+','LineWidth',3,'MarkerSize',8), hold on,
   end


   if 1==0
   figure(10*Mac.plot_Bn2 + 2)
   pcolor(mk,rct,abs(B1PESTN)'), hold on, shading interp
   contour(mk,rct,abs(B1PESTN)','k-'), hold on,
   axis([min(mk) max(mk) 0 1])
   xlabel('m','FontSize',14)
   ylabel('sqrt(\psi_p)','FontSize',14)
   title('|b_s| [Gauss/kA] in PEST CS','FontSize',12)
   colorbar, map = colormap; N = size(map,1); 
   a = linspace(0,1,N)';  b = a.^0.5*ones(1,3); 
   %colormap(map.*b)
   colormap(hot)
   
   figure(10*Mac.plot_Bn2 + 4)
   pcolor(mk,rct,abs(BnGEOMN)'), hold on, shading interp
   contour(mk,rct,abs(BnGEOMN)','k-'), hold on,
   axis([min(mk) max(mk) 0 1])
   xlabel('m','FontSize',14)
   ylabel('sqrt(\psi_p)','FontSize',14)
   title('|b_n| [Gauss/kA] in GEOM. CS','FontSize',12)
   colorbar, map = colormap; N = size(map,1); 
   a = linspace(0,1,N)';  b = a.^0.5*ones(1,3); 
   %colormap(map.*b)
   colormap(hot)
   end

   if CheckPest == 1
   figure(10*Mac.plot_Bn2 + 1)
   pcolor(mk,rct,abs(BnPEST)'), hold on, shading interp
   contour(mk,rct,abs(BnPEST)','k-'), hold on,
   axis([-15 15 0 1])
   xlabel('m','FontSize',14)
   ylabel('sqrt(\psi_p)','FontSize',14)
   title('|b_n| [Gauss/kA] in PEST CS','FontSize',12)
   colorbar, map = colormap; N = size(map,1); 
   a = linspace(0,1,N)';  b = a.^0.5*ones(1,3); 
   %colormap(map.*b)
   colormap(hot)
   
   figure(10*Mac.plot_Bn2 + 3)
   pcolor(mk,rct,abs(B1PEST)'), hold on, shading interp
   contour(mk,rct,abs(B1PEST)','k-'), hold on,
   axis([-15 15 0 1])
   xlabel('m','FontSize',14)
   ylabel('sqrt(\psi_p)','FontSize',14)
   title('|b_s| [Gauss/kA] in PEST CS','FontSize',12)
   colorbar, map = colormap; N = size(map,1); 
   a = linspace(0,1,N)';  b = a.^0.5*ones(1,3); 
   %colormap(map.*b)
   colormap(hot)
   
   figure(10*Mac.plot_Bn2 + 5)
   pcolor(mk,rct,abs(BnGEOM)'), hold on,shading interp
   contour(mk,rct,abs(BnGEOM)','k-'), hold on,
   axis([-15 15 0 1])
   xlabel('m','FontSize',14)
   ylabel('sqrt(\psi_p)','FontSize',14)
   title('|b_n| [Gauss/kA] in GEOM. CS','FontSize',12)
   colorbar, map = colormap; N = size(map,1); 
   a = linspace(0,1,N)';  b = a.^0.5*ones(1,3); 
   %colormap(map.*b)
   colormap(hot)
   end
end


if Mac.plot_Bn3 > 0
   figure(10*Mac.plot_Bn3 + 1)
   plot(rct,real(BnPESTN(1:end,:))), hold on,

   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq],[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 

   figure(10*Mac.plot_Bn3 + 2)
   plot(rct,imag(BnPESTN(1:end,:))), hold on,

   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq],[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 

   figure(10*Mac.plot_Bn3 + 3)
	plot(rct,abs(BnPESTN(1:end,:)),'LineWidth',1),hold on
   xlabel('\psi_p^{1/2}','FontSize',14)
   ylabel('|b_{n=-3}^m| [G/kA]','FontSize',14)
	ha=get(10*Mac.plot_Bn3 + 3,'CurrentAxes'); set(ha,'FontSize',14)

   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq],[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 

   figure(10*Mac.plot_Bn3 + 4)
	plot(mk,max(abs(BnPESTN),[],2),'r-o','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','r'), hold on,
   xlabel('m','FontSize',14)
   ylabel('max|b_{n=-3}^m(\psi_p)| [G/kA]','FontSize',14)
	ha=get(10*Mac.plot_Bn3 + 4,'CurrentAxes'); set(ha,'FontSize',14)
   
end
  

  



function [B1,B2,B3,Bn] = MacGetB123(BM1,BM2,BM3,R,Z,dRdchi,dZdchi)

global Mac SDIR

expmchi = exp(Mac.Mm*Mac.chi*i);

B1 = BM1*expmchi;
B2 = BM2*expmchi;
B3 = BM3*expmchi;

%work out Bn
G22  = dRdchi.^2 + dZdchi.^2;  G22(1,:) = G22(2,:);
N    = size(B1,1);
Bn   = B1./sqrt(G22(1:N,:))./R(1:N,:);
expmchi = exp(-Mac.chi'*Mac.Mm'*i);
BMn = Bn*expmchi*(Mac.chi(2)-Mac.chi(1))/2/pi;

%save Bn along uniform geometric angle at Mac.rs surface
if 1==1
[X,II]=min(abs(Mac.s-Mac.rs));
Bns = Bn(II,:);  %[G]

Rax = R(1,1);
Zax = Z(1,1);
Tgm = atan2(Z(II,:)-Zax,R(II,:)-Rax);
[Tgm,II]=unique(Tgm); Bns = Bns(II);
[Tgm,II]=sort(Tgm);   Bns = Bns(II);
Tnw = linspace(-pi,pi,Mac.Nchi); Tnw=Tnw(:);
Tgm = [Tgm-2*pi Tgm Tgm+2*pi];
Bns = [Bns Bns Bns];
[Tgm,II]=unique(Tgm); Bns=Bns(II);
Bnw = spline(Tgm,Bns,Tnw);
res = [Tnw real(Bnw) imag(Bnw)];
save dBnormal.txt res -ascii -double

%figure(99)
%plot(Tnw,real(Bnw),'r-',Tnw,imag(Bnw),'b-')
end

%plot Bn
if (Mac.plot_Bn>0)
   hf=figure(10*Mac.plot_Bn+0);
   II=1:Mac.Ns1;
   plot(Mac.s(II),real(BMn(II,:)),'LineWidth',1), hold on,
   ylabel('Re[B_n^{(m)}]','FontSize',16,'FontWeight','Bold')
   xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   axis tight
   
   KRFX=0;
   if KRFX==1
       mm = Mac.Mm';
       A   = 4*2; m1=1-mm(1)+1;
       res = A.^abs(mm'-1);
       resa= [mm' res];

	    I1  = Mac.Ns1; 
       res = abs(BMn(I1,m1)./BMn(I1,:))';
       resa= [resa res];
	    [X,I1] = min(abs(Mac.s-1.1)); 
       res = abs(BMn(I1,m1)./BMn(I1,:))';
       resa= [resa res];

       AmpRatioBn= resa
       Bnm1      = BMn(I1,m1)
   end      

   hf=figure(10*Mac.plot_Bn+1);
   C=max(max(abs(real(BMn(II,:))),[],1));
   C1=max(real(BMn(II,:)),[],1);
   C2=min(real(BMn(II,:)),[],1);
   C3=C1;
   I=find(abs(C1)<abs(C2));  C3(I)=C2(I);
   plot(Mac.Mm,abs(C3),'b-o','LineWidth',1,'MarkerSize',8), hold on,
   ylabel('max[Re(B_n^{(m)}(r))]','FontSize',16,'FontWeight','Bold')
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   %axis tight
   
  % contour plot for Bn 
  figure(10*Mac.plot_Bn+2)
  MM = [-3:9] - Mac.Mm(1) + 1;
  pcolor(Mac.Mm(MM),Mac.s(II),abs(BMn(II,MM))), hold on
  shading interp
  contour(Mac.Mm(MM),Mac.s(II),abs(BMn(II,MM)),'k-'), hold on,
  axis([min(Mac.Mm(MM)) max(Mac.Mm(MM)) 0 1])
  xlabel('m','FontSize',16,'FontWeight','Bold')
  ylabel('sqrt(\psi_p)','FontSize',16,'FontWeight','Bold')
  title('|b_n| [Gauss/kA] in PEST CS','FontSize',16,'FontWeight','Bold')
  ha = get(10*Mac.plot_Bn+2,'CurrentAxes');
  set(ha,'FontSize',14,'FontWeight','Bold')
  colorbar, colormap(hot)

   hf=figure(10*Mac.plot_Bn+3);
   %II=1:Mac.Ns1;
   NVG=108; II = 1:Mac.Ns1+NVG;
   pcolor(R(II,:)*Mac.R0EXP,Z(II,:)*Mac.R0EXP,abs(Bn(II,:))), hold on, shading interp
   axis equal
   xlabel('R [m]','FontSize',18,'FontWeight','Bold')
   ylabel('Z [m]','FontSize',18,'FontWeight','Bold')
   title('|{\delta}B_n| [Gauss]','FontSize',14)
   colorbar,  colormap(jet)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
   if 1==1
      %location in [m] of the outboard mid-plane Minov probes in MAST-U (at phi=300 deg)
      Rprobe = 1.84325;	
      Zprobe = -2.29e-08;
      plot(Rprobe,Zprobe,'x','MarkerSize',12,'LineWidth',2)
      II=1:Mac.Ns1+110;
      BnProbe = griddata(R(II,:)*Mac.R0EXP,Z(II,:)*Mac.R0EXP,abs(Bn(II,:)),Rprobe,Zprobe)
   end
   II = Mac.Ns1;
   plot(R(II,:)*Mac.R0EXP,Z(II,:)*Mac.R0EXP,'k-','LineWidth',2), 
   %MacPlotLimiter(10*Mac.plot_Bn+3)

   %get some numbers for Bn
   dBnMaxPls = max(max(abs(Bn(1:Mac.Ns1,:))))
   [Y,J] = max(R(Mac.Ns1,:));
   dBnLFSmid = abs(Bn(Mac.Ns1,J))
end

%plot B1,B2,B3 on the first wall surface along chi
if Mac.plot_Bs > 0
  fac = 1.0;
  [smin,II] = min(abs(Mac.s-Mac.rw(1)));
  figure(Mac.plot_Bs)
  SS = Mac.SS;
  subplot(3,2,1), plot(Mac.chi,real(B1(II,:)),SS,'LineWidth',2), hold on 
                  ylabel('Re[B1]','FontSize',16),
  subplot(3,2,2), plot(Mac.chi,imag(B1(II,:)),SS,'LineWidth',2), hold on 
                  ylabel('Im[B1]','FontSize',16),
  subplot(3,2,3), plot(Mac.chi,real(B2(II,:)),SS,'LineWidth',2), hold on 
                  ylabel('Re[B2]','FontSize',16),
  subplot(3,2,4), plot(Mac.chi,imag(B2(II,:)),SS,'LineWidth',2), hold on 
                  ylabel('Im[B2]','FontSize',16),
  subplot(3,2,5), plot(Mac.chi,real(B3(II,:)),SS,'LineWidth',2), hold on 
                  ylabel('Re[B3]','FontSize',16),
                  xlabel('\chi','FontSize',16),
  subplot(3,2,6), plot(Mac.chi,imag(B3(II,:)),SS,'LineWidth',2), hold on 
                  ylabel('Im[B3]','FontSize',16),
                  xlabel('\chi','FontSize',16),
end


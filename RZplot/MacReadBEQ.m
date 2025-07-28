function MacReadBEQ(NCHIM,R,Z,filename1,filename2)

global Mac

d1 = load(filename1);
d2 = load(filename2);

B1R = d1(:,1);
B1Z = d1(:,2);
B1P = d1(:,3);

B2R = d2(:,1);
B2Z = d2(:,2);
B2P = d2(:,3);

N1S = round(length(B1R)/NCHIM);
N2S = round(length(B2R)/NCHIM);

B1R = reshape(B1R,NCHIM,N1S);
B1Z = reshape(B1Z,NCHIM,N1S);
B1P = reshape(B1P,NCHIM,N1S);

B2R = reshape(B2R,NCHIM,N2S);
B2Z = reshape(B2Z,NCHIM,N2S);
B2P = reshape(B2P,NCHIM,N2S);

BR = transpose([B1R B2R]);
BZ = transpose([B1Z B2Z]);
BP = transpose([B1P B2P]);

BR = [BR BR(:,1)];
BZ = [BZ BZ(:,1)];
BP = [BP BP(:,1)];

NS = size(BR,1);

RR = R(2:NS+1,:)*Mac.R0EXP;
ZZ = Z(2:NS+1,:)*Mac.R0EXP;

Bpol = sqrt(BR.^2+BZ.^2);
BpolMax = max(max(Bpol(1:Mac.Ns1,:)))*Mac.B0EXP;
BpolRes = [Mac.B0EXP BpolMax]

if Mac.plot_BEQ
   hf=figure(10*Mac.plot_BEQ+1);
   pcolor(RR,ZZ,BR*Mac.B0EXP),  colorbar, hold on, shading interp
   plot(R(Mac.Ns1,:)*Mac.R0EXP,Z(Mac.Ns1,:)*Mac.R0EXP,'k-','LineWidth',2), hold on
   xlabel('R [m]','FontSize',16,'FontWeight','Bold')
   ylabel('Z [m]','FontSize',14)
   title('B_R [T]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',12)
   colormap(hot)
   axis equal

   hf=figure(10*Mac.plot_BEQ+2);
   pcolor(RR,ZZ,BZ*Mac.B0EXP),  colorbar, hold on, shading interp
   plot(R(Mac.Ns1,:)*Mac.R0EXP,Z(Mac.Ns1,:)*Mac.R0EXP,'k-','LineWidth',2), hold on
   xlabel('R [m]','FontSize',16,'FontWeight','Bold')
   ylabel('Z [m]','FontSize',14)
   title('B_Z [T]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',12)
   colormap(hot)
   axis equal

   hf=figure(10*Mac.plot_BEQ+3);
   pcolor(RR,ZZ,BP*Mac.B0EXP),  colorbar, hold on, shading interp
   plot(R(Mac.Ns1,:)*Mac.R0EXP,Z(Mac.Ns1,:)*Mac.R0EXP,'k-','LineWidth',2), hold on
   xlabel('R [m]','FontSize',16,'FontWeight','Bold')
   ylabel('Z [m]','FontSize',14)
   title('B_{\phi} [T]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',12)
   colormap(hot)
   axis equal


   hf=figure(10*Mac.plot_BEQ+4);
   BB=sqrt(BR.^2+BZ.^2+BP.^2);   
   LAMBDA=1.1;
   pcolor(RR,ZZ,BB*Mac.B0EXP),  colorbar, hold on, shading interp, colormap(hot)
   plot(R(Mac.Ns1,:)*Mac.R0EXP,Z(Mac.Ns1,:)*Mac.R0EXP,'k-','LineWidth',2), hold on
   [c,h]=contour(RR,ZZ,BB*Mac.B0EXP,[1/LAMBDA 1/LAMBDA]*Mac.B0EXP);
   set(h,'LineWidth',2,'Color','w') 
   xlabel('R [m]','FontSize',16,'FontWeight','Bold')
   ylabel('Z [m]','FontSize',14)
   title('|B| [T]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',12)
   axis equal

   hf=figure(10*Mac.plot_BEQ+5);
   NS = Mac.Ns1-1;
   BB = BB(1:NS,:);
   chin = Mac.chi;
   II = find(chin>pi);
   chin(II) = chin(II)-2*pi;
   [chin,II]=sort(chin);
   BB = BB(:,II);
   [chin,II] = unique(chin);
   BB = BB(:,II);
   chin = chin*180/pi;
   pcolor(chin,Mac.s(2:NS+1).^2,BB*Mac.B0EXP),  colorbar, hold on, shading interp, colormap(hot)
   [c,h]=contour(chin,Mac.s(2:NS+1).^2,BB*Mac.B0EXP,[1/LAMBDA 1/LAMBDA]*Mac.B0EXP);
   %pcolor(chin,Mac.s(2:NS+1).^2,1-BB),  
   %[c,h]=contour(chin,Mac.s(2:NS+1).^2,1-BB,[0 0]);
   if 1==0
   d=load('TEMP_02a');
   d =d(:,[1:4 6:7]);
   II = [9:9]*5;
   s0=d(II+2,:); 
   c0=d(II+3,:); IC=find(c0<0); c0(IC)=c0(IC)+2*pi; 
   plot(c0,s0,'k+',c0(:,5),s0(:,5),'ko','LineWidth',2,'MarkerSize',9)
   end
   set(h,'LineWidth',2,'Color','w') 
   xlabel('poloidal angle [deg.]','FontSize',16,'FontWeight','Bold')
   ylabel('\psi_p','FontSize',14)
   title('|B| [T]','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',12)
   colorbar, shading interp, colormap(hot)

   if 1==0
   hf=figure(10*Mac.plot_BEQ+6);
   tt = d(II+1,5);
   ss = d(II+2,5);
   cc = d(II+3,5);
   yy = d(II+5,5);
   x1 = d(II+1,3);
   x2 = d(II+1,6);
   subplot(3,1,1), plot(tt,ss,'b-+'), ylabel('s'), hold on, a=axis;
                   plot([x1 x1],[a(3) a(4)],'k-','LineWidth',2),
                   plot([x2 x2],[a(3) a(4)],'k-','LineWidth',2)
   subplot(3,1,2), plot(tt,cc,'b-+'), ylabel('\chi'), hold on, a=axis;
                   plot([x1 x1],[a(3) a(4)],'k-','LineWidth',2),
                   plot([x2 x2],[a(3) a(4)],'k-','LineWidth',2)
   subplot(3,1,3), plot(tt,yy,'b-+'), ylabel('1-{\lambda}b'), hold on, a=axis;
                   plot([x1 x1],[a(3) a(4)],'k-','LineWidth',2),
                   plot([x2 x2],[a(3) a(4)],'k-','LineWidth',2)

   end
end

   




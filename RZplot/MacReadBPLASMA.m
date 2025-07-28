function [BM1,BM2,BM3] = MacReadBPLASMA(filename)

global Mac

BPLASMA = load(filename);

Mac.Nm1 = BPLASMA(1,1);

if Mac.Ns > BPLASMA(1,2),
   disp('Number of radial points is different for equilibrium and stability!')
   Ns = BPLASMA(1,2);
   Mac.Ns = Ns;
   Mac.s  = Mac.s(1:Ns);
else 
   Ns = BPLASMA(1,2);
end

Mac.n  = BPLASMA(1,3);
Mac.Mm = round(BPLASMA(2:Mac.Nm1+1,1));

if length(Mac.mm_plot)>0 
   Mac.mm_true = Mac.mm_plot-Mac.Mm(1)+1;
else
   Mac.mm_true = Mac.Mm-Mac.Mm(1)+1;
end

BM1 = BPLASMA(Mac.Nm1+2:end,1) + BPLASMA(Mac.Nm1+2:end,2)*i;
BM2 = BPLASMA(Mac.Nm1+2:end,3) + BPLASMA(Mac.Nm1+2:end,4)*i;
BM3 = BPLASMA(Mac.Nm1+2:end,5) + BPLASMA(Mac.Nm1+2:end,6)*i;

BM1 = reshape(BM1,Ns,Mac.Nm1);
BM2 = reshape(BM2,Ns,Mac.Nm1);
BM3 = reshape(BM3,Ns,Mac.Nm1);

BM1 = BM1(1:Mac.Ns,:);
BM2 = BM2(1:Mac.Ns,:);
BM3 = BM3(1:Mac.Ns,:);

% output BM1 at rw-surface
if 1==1
[smin,II] = min(abs(Mac.s-Mac.rw(1)));
res = BM1(II,:);
res = res(:);
res = [Mac.Nm1 II; real(res) imag(res)];
save INPUT_BNM res -ascii -double
end

BM1 = BM1*Mac.BNORM;  %[T]
BM2 = BM2*Mac.BNORM;  %[T]
BM3 = BM3*Mac.BNORM;  %[T]

%patch BM2 near magnetic axis
if 1==1
   BM2(1,:) = 0;
   BM2(2,:) = BM2(3,:)/2;
end

if 1==0
   N = Mac.Ns1;
   tmp = real(BM1(1:N,:));
   save b1m_real.txt tmp -ascii 
   tmp = imag(BM1(1:N,:));
   save b1m_imag.txt tmp -ascii 
   tmp = Mac.s(1:N).^2;
   save psip.txt tmp -ascii
end 

% subtract combined basis solution
% for Delta' computation
if Mac.FullSol==1
   load temp_B1ATT
   BM1 = -BM1;
   BM1(Mac.Iratsurf,:) = transpose(B1ATT) + BM1(Mac.Iratsurf,:);
   res_BM1 = [Mac.Mm transpose(BM1(Mac.Iratsurf,:))]
end


% convergence extraplation
global knorm
knorm = 0;
if knorm == 1
[smin,II] = min(abs(Mac.s-Mac.rw(2)));
BM10 = sum(BM1(II+1:end-1,:),2);
for m=1:Mac.Nm1
  BM1(II+1:end-1,m) = BM1(II+1:end-1,m)./BM10;
  BM1(II,m) = spline(Mac.s(II+1:end-1),BM1(II+1:end-1,m),Mac.s(II));
end
end

% test: for Delta' calculations
if 1 == 0
   data = load('/.automount/funsrv1/root/home/yliu/DeltapAna/PROFEQ.OUT');
   q = data(:,2);
   II = 1:Mac.Ns1;
   for k=1:size(BM1,2)
       BM1(II,k) = BM1(II,k)./q;
   end
end

% save BNORM for MARS-F input file, used in BOWALL02 or BOWALL03
if (length(Mac.rw) > 0)
  [smin,II] = min(abs(Mac.s-Mac.rw(1)));
  if (Mac.save_BNORM==1)
    tmp = BM1(II,:)/Mac.BNORM;
    BNORM = [real(tmp(:)) imag(tmp(:))];
    save BNORM01_MARSF BNORM -ascii -double
  elseif (Mac.save_BNORM==2)
    tmp = BM1(II,:);
    BNORM = tmp/BM1(II,1);
    BNORM = [real(BNORM(:)) imag(BNORM(:))];
    save BNORM02_MARSF BNORM -ascii -double
  elseif (Mac.save_BNORM==3)
    tmp = BM1(II,:);
    BNORM = [real(tmp(:)) imag(tmp(:))];
    tmp = BM2(II-1,:);
    BNORM = [BNORM real(tmp(:)) imag(tmp(:))];
    tmp = BM1(II,:)./BM2(II-1,:);
    BNORM = [real(tmp(:)) imag(tmp(:))];
    save BNORM03_MARSF BNORM -ascii -double
  elseif (Mac.save_BNORM==4)
    tmp = BM1(II,:)/sum(BM2(II-1,:));
    BNORM = [real(tmp(:)) imag(tmp(:))];
    save BNORM04_MARSF BNORM -ascii -double
  end
end

% note that BM2 and BM3 are defined at half-points, recompute at integer-points
if Mac.spline_B23 == 1
  x = (Mac.s(1:Mac.Ns-1) + Mac.s(2:Mac.Ns))*0.5;
  BM2new = BM2;
  BM2new(2:end-1,:) = transpose(pchip(x',transpose(BM2(1:end-1,:)),Mac.s(2:Mac.Ns-1)'));
  BM2new(1,:) = 0;  BM2new(end,:) = BM2new(end-1,:);  BM2 = BM2new;
  BM3new = BM3;
  BM3new(2:end-1,:) = transpose(pchip(x',transpose(BM3(1:end-1,:)),Mac.s(2:Mac.Ns-1)'));
  BM3new(1,:) = 0;  BM3new(end,:) = BM3new(end-1,:);  BM3 = BM3new;
  clear BM2new BM3new
elseif Mac.spline_B23 == 2
  BM2(2:end,:) = BM2(1:end-1,:);
  BM3(2:end,:) = BM3(1:end-1,:);
elseif Mac.spline_B23 == 3
  BM2(end,:) = transpose(pchip(Mac.s(1:end-1)',transpose(BM2(1:end-1,:)),Mac.s(end)));
  BM3(end,:) = transpose(pchip(Mac.s(1:end-1)',transpose(BM3(1:end-1,:)),Mac.s(end)));
end

if Mac.plot_BM > 0
   figure(10*Mac.plot_BM+0)
   SS = '-';
   mm = Mac.mm_plot; 
   JJ=mm-Mac.Mm(1)+1;
   if JJ(end)>length(Mac.Mm), JJ=[1:length(Mac.Mm)]; end   
   BM = BM1(:,JJ); a = angle(BM); C = max(max(abs(BM)));
   I = find(abs(BM)/C<0.2); a(I)=0; I = find(a<0); a(I)=a(I)+pi;
   %subplot(3,2,1), plot(Mac.s,a/pi,SS), hold on,
   subplot(3,2,1), plot(Mac.s,real(BM),SS), hold on,
                   ylabel('Re(B^1_m)','FontSize',14)
   subplot(3,2,2), plot(Mac.s,imag(BM),SS), hold on,
                   xlabel('s','FontSize',14), ylabel('Im(B^1_m)','FontSize',14)

   BM = BM2(:,JJ); a = angle(BM); C = max(max(abs(BM)));
   I = find(abs(BM)/C<0.2); a(I)=0; I = find(a<0); a(I)=a(I)+pi;
   subplot(3,2,3), plot(Mac.s,real(BM),SS), hold on,
                   ylabel('Re(B^2_m)','FontSize',14)
   subplot(3,2,4), plot(Mac.s,imag(BM),SS), hold on,
                   xlabel('s','FontSize',14), ylabel('Im(B^2_m)','FontSize',14)

   BM = BM3(:,JJ); a = angle(BM); C = max(max(abs(BM)));
   I = find(abs(BM)/C<0.2); a(I)=0; I = find(a<0); a(I)=a(I)+pi;
   subplot(3,2,5), plot(Mac.s,real(BM),SS), hold on,
                   ylabel('Re(B^3_m)_m','FontSize',14)
   subplot(3,2,6), plot(Mac.s,imag(BM),SS), hold on,
                   xlabel('s','FontSize',14), ylabel('Im(B^3_m)','FontSize',14)

   KDELTAP = 0;

   hff=figure(10*Mac.plot_BM+1);
   NsII = min(Mac.Ns,Mac.Ns1+10);
   II=[1:Mac.Ns1];
   %mm=Mac.Mm;
   mm=Mac.mm_plot;
   JS = 50;
   if length(Mac.Iratsurf)>0
      JJ=Mac.Iratsurf;
      JS=JJ(Mac.Kratsurf);
   end

   if KDELTAP==2, II=[Mac.Iratsurf+1:Mac.Ns1]; end      
   if KDELTAP==-1, II=[1:Mac.Ns1+60]; mm=[4 5]; end
   if KDELTAP==3, II=[JS+1:JS+10]; end      
   if KDELTAP==4, II=[JS+1:JS+20]; end      
   if KDELTAP==5, II=[JS+1:JS+20]; end      
   if KDELTAP==5, II=[JS+1:JS+20]; end      
   if KDELTAP==6, NS=15; II=[[JS-NS:JS-1] [JS+1:JS+NS]]; end      

   JJ=mm-Mac.Mm(1)+1;
   if JJ(end)>length(Mac.Mm), JJ=[1:length(Mac.Mm)]; end   

   Y = abs(BM1(II,JJ)); %Cn=max(max(abs(Y))); Y = Y/Cn;


   if KDELTAP==1
   xo = Mac.s(II); yo = Y;
   xx = (xo(1:end-1)+xo(2:end))/2; 
   M = size(yo,2); yx=zeros(length(xx),M);
   for k=1:M
       yx(:,k) = diff(yo(:,k))./diff(xo);
   end
   hp=plot(xx,yx,'-','LineWidth',2); hold on,
   ylabel('Re[\partial{b^1}/\partial{s}]','FontSize',16,'FontWeight','Bold')
   end

   if KDELTAP==11
   nu0 = 9.0383e-03;
   nu0 = 0;
   xo = Mac.s(II); yo = Y;
   xs = Mac.s(JS);
   M = size(yo,2); yx=yo*0;
   IK = find(xo ~= xs);
   for k=1:M
       yx(IK,k) = yo(IK,k)./(abs(xo(IK)-xs)).^(1+nu0);
   end
   hp=plot(xo,yx,'-','LineWidth',2); hold on,
   ylabel('Re[\partial{b^1}/\partial{s}]','FontSize',16,'FontWeight','Bold')
   end

   if KDELTAP==0 | KDELTAP==-1
   hp=plot(Mac.s(II),Y,'k-','LineWidth',3); hold on,
   if length(Mac.Iratsurf)>0
      sres = Mac.s(Mac.Iratsurf(1));
      a=axis; plot([sres sres],[a(3) a(4)],'k--')
   end
   ylabel('|b^1=qb\cdot\nabla\psi_p/(B\cdot\nabla\phi)|','FontSize',16,'FontWeight','Bold')
   res = [Mac.s(II) Y];
   save FAR_BM.txt res -ascii    
   end
   if KDELTAP==2
   xo = Mac.s(II)-Mac.s(JS); yx = Y;
   yo = zeros(size(yx,1)-1,size(yx,2));
   for k=1:size(yx,2)
       yo(:,k) = diff(yx(:,k))./diff(xo);
   end
   xo = (xo(1:end-1)+xo(2:end))/2;
   hp=loglog(xo,abs(yo),'-','LineWidth',2); hold on,
   ylabel('|Re[b^1]|','FontSize',16,'FontWeight','Bold')

   nu = (log(abs(yo(2,:)))-log(abs(yo(1,:))))./(log(xo(2))-log(xo(1)));
   logDelta = log(abs(yo(1,:))) - (0+nu)*log(xo(1));
   Delta = exp(logDelta);
   IK = find(yo(1,:)<0); Delta(IK)=-Delta(IK);
   yn = yo;
   for k=1:size(yo,2)
       yn(:,k) = abs(Delta(k))*xo.^(0+nu(k));
   end
   loglog(xo,yn,'--'), hold on,
   Delta = Delta./(1+nu);
   res_Delta = [nu' Delta']
   end
   if KDELTAP==3
      hp=plot(Mac.s(II),Y,'-','LineWidth',2); hold on,
      ylabel('Re[b^1=qb\cdot\nabla\psi_p/(B\cdot\nabla\phi)]','FontSize',16,'FontWeight','Bold')
      xo = Mac.s(II)-Mac.s(JS); yo = Y;
      [nu,Delta,Fmin] = MacDeltapGetNu(xo,yo,2);
      res_Delta = [nu Delta Fmin]
   end
   if KDELTAP==4
      xs = 6.73106839e-1;
      xo = Mac.s(II)-xs; yo = Y;
      xx = (xo(1:end-1)+xo(2:end))/2; 
      M = size(yo,2); yx=zeros(length(xx),M);
      for k=1:M
          yx(:,k) = diff(yo(:,k))./diff(xo);
      end
      hp=plot([xs;Mac.s(II)],[zeros(1,size(Y,2)); Y],'-+','LineWidth',2); hold on,
      ylabel('Re[b^1=qb\cdot\nabla\psi_p/(B\cdot\nabla\phi)]','FontSize',16,'FontWeight','Bold')
      nu0 = 2.6324e-1;
      Delta = MacDeltapGetD(xx,yx,nu0,1,[2 4],1);
      %Delta = MacDeltapGetDls(xx,yx,nu0,[1:5],1);
      res_Delta = Delta/(1+nu0)
   end
   if KDELTAP==5
      xs = Mac.s(JS);
      xo = Mac.s(II)-xs; yo = Y;
      hp=plot([xs;Mac.s(II)],[zeros(1,size(Y,2)); Y],'-+','LineWidth',2); hold on,
      ylabel('Re[b^1=qb\cdot\nabla\psi_p/(B\cdot\nabla\phi)]','FontSize',16,'FontWeight','Bold')
      nu0 = 1.3381e-02;
      Delta = MacDeltapGetDs(xo,real(yo),nu0,2,[1:8],1);
      res_Delta = Delta
   end
   if KDELTAP==6
      xs = Mac.s(JS);
      xo = Mac.s(II)-xs; yo = Y;
      hp=plot(Mac.s(II),Y,'-+','LineWidth',2); hold on,
      ylabel('Re[b^1=qb\cdot\nabla\psi_p/(B\cdot\nabla\phi)]','FontSize',16,'FontWeight','Bold')
      nu0 = 3.1187e-02;
      %Delta = MacDeltapGetDf(xo,yo,nu0,[NS-13:NS+14],1);
      Delta = MacDeltapGetDff(xo,yo,nu0,[NS-13:NS+14],1);
      res_Delta = Delta
   end

   figure(hff)
   xlabel('s\equiv\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hff,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   for k=1:length(mm)
       c = get(hp(k),'Color');
       [X,I]=max(abs(Y(:,k)));
       text(Mac.s(I),Y(I,k),int2str(mm(k)),'FontSize',18,'FontWeight','Bold','Color',c)
   end
   %text(0.05,0.95,'(c)','FontSize',18,'FontWeight','Bold')
   axis tight

   KRFX=0;
   if KRFX==1
       A   = 4*2; m1=1-mm(1)+1;
       res = A.^abs(mm'-1);
       resa= [mm' res];

	    II  = Mac.Ns1; 
       res = real(BM1(II,m1)./BM1(II,:))';
       resa= [resa res];
	    [X,II] = min(abs(Mac.s-1.1)); 
       res = real(BM1(II,m1)./BM1(II,:))';
       resa= [resa res];

       AmpRatioB1= resa;
       B1m1      = BM1(II,m1)
   end      

   hf=figure(10*Mac.plot_BM+2);
   C=max(max(abs(abs(BM1(1:Mac.Ns1,:))),[],1));
   C1=max(abs(BM1(1:Mac.Ns1,:)),[],1);
   C2=min(abs(BM1(1:Mac.Ns1,:)),[],1);
   C3=C1;
   I=find(abs(C1)<abs(C2));  C3(I)=C2(I);
   plot(Mac.Mm,C3,'b-o','LineWidth',1,'MarkerSize',8), hold on,
   ylabel('max[Abs(b^1_{(m)}(r))]','FontSize',16,'FontWeight','Bold')
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   axis tight
end


% go through database 
% plot SVD-MoR data

%SDIR1 = '/cscratch/liuy/WorkIR_D/Database/';
SDIR1 = '/home/liuy/Work/IR_D/Database/';
SDIRM = '/home/liuy/IRIS/mars-code/RZplot/';

%DEV = {'MAST','DIII-D','AUG','ITER'};
%DEV = {'D3D'};
DEV = {'DIII-D'};
SN  = '1234';
LN  = zeros(1,size(DEV,2));
SC  = 'U';
SM  = {'g+', 'ro', 'bs', 'kd'};
ST  = [50 20 10];
kd1 = 1;
kSUV = 0;  %plot first kSUV components 

S1J = []; 
S2J = []; 
SNJ = []; 
REJ = [];
REJR= [];
S1X = []; 
S2X = []; 
SNX = []; 
REX = [];
REXR= [];
INDX= [];

for kd=kd1:size(DEV,2)
    if kd==1 | kd==3, SC='ULC';  end
    %if kd==1 | kd==3, SC='U';  end
    if kd==2,         SC='ULC'; end
    if kd==4,         SC='ULM'; end

    SDIR2 = [SDIR1 DEV{kd} '/'];
    cd(SDIR2)
    F  = split(ls);
    NF = length(F)-1
    cd(SDIRM)

    for kf=1:NF
    %for kf=23:23
        FF = F{kf};
        %for kn=1:length(SN)
        for kn=4:4
        for kc=1:length(SC)
        %for kc=1:1
            SDIR3 = [SDIR2 FF '/3D/n' SN(kn) '/'];
            SDIR4 = [SDIR3 SC(kc) '/'];
            dj    = load([SDIR4 'SVD_Jpara5.txt']);
            dx    = load([SDIR4 'SVD_X1U5.txt']);

            N   = round(size(dj,2)/2);
            S1J = [S1J; dj(2,1)];
            S2J = [S2J; dj(2,2)];
            SNJ = [SNJ; dj(2,N)];
            REJR= [REJR; dj(2,end-1)];
            REJ = [REJ; dj(2,end)];

            S1X = [S1X; dx(2,1)];
            S2X = [S2X; dx(2,2)];
            SNX = [SNX; dx(2,N)];
            REXR= [REXR; dx(2,end-1)];
            REX = [REX; dx(2,end)];

            INDX = [INDX; kd];

            if kSUV>0
               dp = load([SDIR3 'PROFEQ.OUT']);
               s  = dp(:,1); 
               q  = dp(:,2);
               mmin = ceil(3*min(q));
               mmax = floor(3*max(q));
               NS = length(s);
               m  = [-33:33];
               NM = length(m);
 
               hf = figure(11);
               semilogy(dj(2,1:kSUV),'+'), hold on,
               SVD_S = dj(2,1:kSUV)

               hf = figure(12);
               II = 3:NS+2;
               JJ = 1:kSUV;
               XL=0.2; XW=0.75; YB=0.15; YW=0.8; AW=0.1; N=2;
               hs=subplot('position',[XL YB+YW*(N-1)/N XW YW/N]); 
               plot(s.^2,dj(II,JJ),'-','LineWidth',2), hold on
               
               hf = figure(12);
               II = 3:NS+2;
               JJ = N + [1:kSUV];
               hs=subplot('position',[XL YB+YW*(N-2)/N XW YW/N]); 
               plot(s.^2,dj(II,JJ),'-','LineWidth',2), hold on
               
               hf = figure(14);
               II = NS+3:NS+3+NM-1;
               JJ = 1:kSUV;
               hs=subplot('position',[XL YB+YW*(N-1)/N XW YW/N]); 
               plot(m,dj(II,JJ),'-','LineWidth',2), hold on
               
               hf = figure(14);
               II = NS+3:NS+3+NM-1;
               JJ = N + [1:kSUV];
               hs=subplot('position',[XL YB+YW*(N-2)/N XW YW/N]); 
               plot(m,-dj(II,JJ),'-','LineWidth',2), hold on

               hf = figure(21);
               semilogy(dx(2,1:kSUV),'+'), hold on,

               hf = figure(22);
               II = 3:NS+2;
               JJ = 1:kSUV;
               plot(s,dx(II,JJ)), hold on
               
               hf = figure(23);
               II = 3:NS+2;
               JJ = N + [1:kSUV];
               plot(s,dx(II,JJ)), hold on
               
               hf = figure(24);
               II = NS+3:NS+3+NM-1;
               JJ = 1:kSUV;
               plot(m,dx(II,JJ)), hold on
               
               hf = figure(25);
               II = NS+3:NS+3+NM-1;
               JJ = N + [1:kSUV];
               plot(m,dx(II,JJ)), hold on
            end
         end
         end
     end
end

if kSUV > 0
   hf=figure(11);
   xlabel(['SVD mode number'],'FontSize',18,'FontWeight','Bold')
   ylabel('S-value','FontSize',18,'FontWeight','Bold')
   title('{\delta}J_{||}','FontSize',14,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

   hf=figure(12);
   XL=0.2; XW=0.75; YB=0.15; YW=0.8; AW=0.1; N=2;
   hs=subplot('position',[XL YB+YW*(N-1)/N XW YW/N]); 
   %xlabel(['\psi_p'],'FontSize',18,'FontWeight','Bold')
   ylabel('Re(U)','FontSize',18,'FontWeight','Bold')
   %title('{\delta}J_{||}','FontSize',14,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
   %legend('Re[U(1)]','Re[U(2)]','Re[U(3)]','Re[U(4)]','Re[U(5)]')
   a=axis; axis([0.8 1 a(3) a(4)])
   set(hs,'XTick',[])
   for m=mmin:mmax
       [Y,II]=min(abs(q-m/3));
       plot([s(II)^2 s(II).^2],[a(3) a(4)],'k--')
   end

   hf=figure(12);
   hs=subplot('position',[XL YB+YW*(N-2)/N XW YW/N]); 
   xlabel(['\psi_p'],'FontSize',18,'FontWeight','Bold')
   ylabel('Im(U)','FontSize',18,'FontWeight','Bold')
   %title('{\delta}J_{||}','FontSize',14,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
   a=axis; axis([0.8 1 a(3) a(4)])
   for m=mmin:mmax
       [Y,II]=min(abs(q-m/3));
       plot([s(II)^2 s(II).^2],[a(3) a(4)],'k--')
   end
   legend('U(1)','U(2)','U(3)','U(4)','U(5)')

   hf=figure(14);
   hs=subplot('position',[XL YB+YW*(N-1)/N XW YW/N]); 
   %xlabel(['m'],'FontSize',18,'FontWeight','Bold')
   ylabel('Re(V)','FontSize',18,'FontWeight','Bold')
   %title('{\delta}J_{||}','FontSize',14,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
   %legend('Re[V(1)]','Re[V(2)]','Re[V(3)]','Re[V(4)]','Re[V(5)]')
   a=axis; axis([-5 33 a(3) a(4)])
   set(hs,'XTick',[])

   hf=figure(14);
   hs=subplot('position',[XL YB+YW*(N-2)/N XW YW/N]); 
   xlabel(['m'],'FontSize',18,'FontWeight','Bold')
   ylabel('Im(V)','FontSize',18,'FontWeight','Bold')
   %title('{\delta}J_{||}','FontSize',14,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
   legend('V(1)','V(2)','V(3)','V(4)','V(5)')
   a=axis; axis([-5 33 a(3) a(4)])
   set(hs,'XTick',[-5 0 10 20 30])

   hf=figure(21);
   xlabel(['SVD mode number'],'FontSize',18,'FontWeight','Bold')
   ylabel('S-value','FontSize',18,'FontWeight','Bold')
   title('\xi^1','FontSize',14,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

   hf=figure(22);
   xlabel(['s'],'FontSize',18,'FontWeight','Bold')
   ylabel('Re(U)','FontSize',18,'FontWeight','Bold')
   title('\xi^1','FontSize',14,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

   hf=figure(23);
   xlabel(['s'],'FontSize',18,'FontWeight','Bold')
   ylabel('Im(U)','FontSize',18,'FontWeight','Bold')
   title('\xi^1','FontSize',14,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

   hf=figure(24);
   xlabel(['m'],'FontSize',18,'FontWeight','Bold')
   ylabel('Re(V)','FontSize',18,'FontWeight','Bold')
   title('\xi^1','FontSize',14,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

   hf=figure(25);
   xlabel(['m'],'FontSize',18,'FontWeight','Bold')
   ylabel('Im(V)','FontSize',18,'FontWeight','Bold')
   title('\xi^1','FontSize',14,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
end

hf=figure(1);
X=SNJ./S1J; Y=REJ*100; 
for kd=kd1:size(DEV,2)
    I=find(INDX==kd); LN(kd)=length(I);
    plot(X(I),Y(I),SM{kd},'MarkerSize',8), hold on
end
xlabel(['S_' int2str(N) '/S_1'],'FontSize',18,'FontWeight','Bold')
ylabel('Relative Error (%)','FontSize',18,'FontWeight','Bold')
title('{\delta}J_{||}','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
axis([0 1 0 100]);
for kt=1:length(ST)
    I=find(Y<=ST(kt)); 
    Z=length(I)/length(Y)*100;
    plot([0 0.8],[ST(kt) ST(kt)],'k--')
    text(0.8,ST(kt),[num2str(Z) '%'],'FontSize',18,'FontWeight','Bold')
end
%legend(DEV{kd1:4})

RES_LN = LN  %=[280   336   112   384]

hf=figure(2);
X=SNX./S1X; Y=REX*100; 
for kd=kd1:size(DEV,2)
    I=find(INDX==kd); 
    plot(X(I),Y(I),SM{kd},'MarkerSize',8), hold on
end
xlabel(['S_' int2str(N) '/S_1'],'FontSize',18,'FontWeight','Bold')
ylabel('Relative Error (%)','FontSize',18,'FontWeight','Bold')
title('{\xi}^1','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
axis([0 1 0 100]);
for kt=1:length(ST)
    I=find(Y<=ST(kt)); 
    Z=length(I)/length(Y)*100;
    plot([0 0.8],[ST(kt) ST(kt)],'k--')
    text(0.8,ST(kt),[num2str(Z) '%'],'FontSize',18,'FontWeight','Bold')
end
%legend(DEV{kd1:4})



% plot boundary and limietr shapes from g-file
% from LIMITER.OUT produced by MARS-F (REORBIT module)

function MacPlotLimiter(kfig)

global SDIR

%SDIR = '/cscratch/liuy/WorkReKink/';
%SDIR = '/cscratch/liuy/WorkD3D177038/';
%SDIR = '/home/liuy/Work/D3D177038/Data_RE/Common/';
%SDIR = '/home/liuy/Work/ITER/Data_RE/Common/';
%SDIR = '/home/liuy/Work/COMPASS/Data/Common/';
%SDIR = '/home/liuy/Work/RE_KINK/Data_REORBIT/';
%SDIR = '/home/liuy/Work/RE_KINK/Data_REK/Common/';
%SDIR = '/home/liuy/Work/JET/Data/';
%SDIR = '~/Work/MAST/DataN/';

%SDIR = '/home/liuy/Work/RE_KINETIC/Data/';

d = load([SDIR 'LIMITER.OUT']);

RLFT=d(1,1); RLEN=d(1,2); RRGT=RLFT+RLEN;
ZLFT=d(2,1); ZLEN=d(2,2); ZRGT=ZLFT+ZLEN;

NB=d(3,1); NL=d(3,2);
RB=d(4:3+NB,1); 
ZB=d(4:3+NB,2); 

RL=d(4+NB:end,1); 
ZL=d(4+NB:end,2);

%RL = RL(20:90);
%ZL = ZL(20:90);

hf = figure(kfig);
plot([RLFT RRGT RRGT RLFT RLFT],[ZLFT ZLFT ZRGT ZRGT ZLFT],'k-','LineWidth',2), hold on,
%plot(RB,ZB,'g--','LineWidth',0.5), hold on,
plot(RL,ZL,'g-','LineWidth',3), hold on,
axis equal
xlabel('R [m]','FontSize',18,'FontWeight','Bold')
ylabel('Z [m]','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

% check whether a given point is inside or outside limiter boundary
if 1==0
   %compute geometric angle, with origin point O=center of plasma boundary
   R0=(max(RB)+min(RB))/2;
   Z0=(max(ZB)+min(ZB))/2;
   
   AL=atan2(ZL-Z0,RL-R0);
   II=find(AL<0); AL(II)=AL(II)+2*pi;

   %for each segment in (RL,ZL), find point X on the segment, 
   %such that OX is perpendicular to the segment
   RX = RL(1:end-1);  ZX = RX;
   for J=1:length(RL)-1
       R1=RL(J);   Z1=ZL(J);
       R2=RL(J+1); Z2=ZL(J+1);
       DEL = (R2-R1)^2+(Z2-Z1)^2;
       R3  = (R2-R1)^2*R0+(Z2-Z1)^2*R1-(R2-R1)*(Z2-Z1)*(Z1-Z0);
       Z3  = (R2-R1)^2*Z1+(Z2-Z1)^2*Z0-(R2-R1)*(Z2-Z1)*(R1-R0);
       RX(J) = R3/DEL;
       ZX(J) = Z3/DEL;
   end

   %compute geometric angle of point X, using the same origin point O
   AX=atan2(ZX-Z0,RX-R0);
   II=find(AX<0); AX(II)=AX(II)+2*pi;

   %give an arbitrary point Y
   RY = 2.5;
   ZY = 0;

   plot(RY,ZY,'rx','LineWidth',1,'MarkerSize',9)

   %compute geometric angle of OY
   AY = atan2(ZY-Z0,RY-R0);
   if (AY<0), AY=AY+2*pi; end

   %find which segment of (RL,ZL) will intersect a line aligned with OY
   II = find((AY-AL(1:end-1)).*(AY-AL(2:end))<0 & (abs(AY-AL(1:end-1))+abs(AY-AL(2:end)))<1.9*pi); 

   INTERSECT=0;
   if (length(II)>0), 
      II=II(1);
      LY = sqrt((RY-R0)^2+(ZY-Z0)^2);
      LX = sqrt((RX(II)-R0)^2+(ZX(II)-Z0)^2); 
      LL = LX/cos(AX(II)-AY);
      if (LY>=LL), INTERSECT=1; end

      plot([RL(II) RL(II+1)],[ZL(II) ZL(II+1)],'ro','MarkerSize',9)
      plot([R0 RX(II)],[Z0 ZX(II)],'r-')
      RY2 = R0 + (RY-R0)*10;
      ZY2 = Z0 + (ZY-Z0)*10;
      a = axis; plot([R0 RY2],[Z0 ZY2],'r--'), axis(a);
   end

   disp(['INTERSECT=',int2str(INTERSECT)])
end
 

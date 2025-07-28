%function MacRfaCtJ3
%% The 2nd (fast) method to
%% transform J3U from MARS coordinates to PEST coordinates
%% and plot thevpoloidal harmonics in PEST
%% before running this subroutine, make sure that 
%% MacMain* is a procedure, not a function
%%   1) RMZM_F & JPLASMA are used in the MacMain* subroutine
%%   2) SDIR defined here agrees with those in the MacMain* subroutine 
%%   3) files RMZM_F_EQAC, RMZM_F_PEST,JXPLASMA(EQAC) exist in SDIR
%%   4) comment out Mac.RunJ in MacMain*

global Mac 
global SDIR 

SDIR = '/home/liuy/IRIS/WorkRWM_TTT/'; SMAIN='TTT_J';

facn      = 1.0;
Mac.RunV  = 0;
Mac.RunB  = 0;
n         = 1;
CheckQ    = 1;
PlotJ3    = 6;
PlotQ     = 5;
LSS       = '-o'; 
%LSC       = [(kclr-1)/(Nclr-1) (kclr-1)*(Nclr-kclr)*4/Nclr^2 (Nclr-kclr)/(Nclr-1)]
LSC       = [0 0 1];

mk = -9:33;
mi = mk - mk(1) + 1;
  
eval(['copyfile ' SDIR 'RMZM_F_EQAC ' SDIR 'RMZM_F.OUT']);
Mac.RunJ = 1;
eval(['MacMain' SMAIN]);
II = 1:Mac.Ns1;
JphiEQAC = Jphi(II,:);
R_EQAC = R(II,:);
Z_EQAC = Z(II,:);

eval(['copyfile ' SDIR 'RMZM_F_PEST ' SDIR 'RMZM_F.OUT']);
Mac.RunJ = 0;
eval(['MacMain' SMAIN]);
ss     = Mac.s(II);
R_PEST = R(II,:);
Z_PEST = Z(II,:);

ss        = Mac.s(II);
JphiPEST  = griddata(R_EQAC,Z_EQAC,JphiEQAC,R_PEST,Z_PEST);
J3PEST    = JphiPEST.*jacobian(II,:)./R_PEST;

expmchi = exp(-Mac.chi'*mk*i);
JM3PEST = J3PEST*expmchi*(Mac.chi(2)-Mac.chi(1))/2/pi;

%patch possible NaN's of JM3PEST near the plasma boundary
if 1==1
for k=1:size(JM3PEST,2)
    INOR = find(isnan(JM3PEST(:,k))==0);
    if length(INOR) < size(JM3PEST,1)
       JM3PEST(:,k) = pchip(ss(INOR),JM3PEST(INOR,k),ss);
    end
end
end
JM3PEST(1,:) = JM3PEST(2,:);

if CheckQ > 0
   dataq = load([SDIR 'PROFEQ.OUT']);
   s = dataq(:,1);
   q = dataq(:,2);
   dpsi = dataq(:,12);
end


if PlotJ3 > 0
   %plot rational surfaces
   if CheckQ > 0
      mq = [ceil(min(q)*abs(n)):max(mk)];
      qq = mq/abs(n);
      [sq,qn] = MacFindX(s,q,qq);
      mq = qn*n;
      mm2 = mq; II2=mm2-mk(1)+1;
   end

   if CheckQ > 0
      hf=figure(10*PlotQ + 0);
      plot(s.^2,q,LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
      xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
      ylabel('q','FontSize',16,'FontWeight','Bold')
      ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
      a = axis;
      for k=1:length(sq)
          plot([0 1],[mq(k) mq(k)]/n,'k--'), hold on,
          plot([sq(k)^2 sq(k)^2],[a(3) a(4)],'k--'), hold on,
      end
   end
   
   %plot real part of resonant harmonics of J3
   hf=figure(10*PlotJ3 + 1);
   plot(ss.^2,real(JM3PEST(:,II2)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('Re(J^3_m) [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq].^2,[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 

   %plot imaginary part of resonant harmonics of J3
   hf=figure(10*PlotJ3 + 2);
   plot(ss.^2,imag(JM3PEST(:,II2)),LSS(1:end-1),'LineWidth',2,'Color',LSC), hold on,
   xlabel('\psi_p','FontSize',16,'FontWeight','Bold')
   ylabel('Im(J^3_m) [mm]','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
   if CheckQ > 0
      a = axis;
      plot([0 1],[0 0],'k-'), hold on,
      plot([sq; sq].^2,[a(3)*ones(size(sq)); a(4)*ones(size(sq))],'k--'), hold on,
   end 

end

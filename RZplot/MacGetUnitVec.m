function [dRds,dZds,dRdchi,dZdchi,jacobian] = MacGetUnitVec(R,Z)

global Mac

dRds = R; dZds = R; dRdchi = R; dZdchi = R; jacobian = R;

%vacuum region
II = Mac.Ns1:size(R,1);
II2 = Mac.Ns1+1:size(R,1);

s0 = transpose(Mac.s(II));   R0 = R(II,:);
chi0 = Mac.chi;          Z0 = Z(II,:);
hs = min(s0(2:end)-s0(1:end-1))/2; hs = min(hs,2e-5);
hchi = min(chi0(2:end)-chi0(1:end-1))/2; hchi = min(hchi,1e-4);
s1 = s0 - hs;        s2 = s0 + hs;  
chi1 = chi0 - hchi;  chi2 = chi0 + hchi;  

% compute dR/ds using R(s,chi) and spline
R1 = spline(s0,R0',s1);  R2 = spline(s0,R0',s2);
dRds(II,:) = transpose((R2-R1)/hs/2); 

% compute dZ/ds using Z(s,chi) and spline
Z1 = spline(s0,Z0',s1);  Z2 = spline(s0,Z0',s2);
dZds(II,:) = transpose((Z2-Z1)/hs/2); 

% compute dR/dchi using R(s,chi) and spline
R1 = spline([chi0(end-1)-2*pi chi0 chi0(2)+2*pi],[R0(:,end-1) R0 R0(:,2)],chi1);  
R2 = spline([chi0(end-1)-2*pi chi0 chi0(2)+2*pi],[R0(:,end-1) R0 R0(:,2)],chi2);
R1(:,1) = R1(:,end);  R2(:,end) = R2(:,1); 
dRdchi(II,:) = (R2-R1)/hchi/2; 


% compute dZ/dchi using Z(s,chi) and spline
Z1 = spline(chi0,Z0,chi1);  Z2 = spline(chi0,Z0,chi2);
Z1(:,1) = Z1(:,end);  Z2(:,end) = Z2(:,1); 
dZdchi(II,:) = (Z2-Z1)/hchi/2; 

%plasma region
II = 1:Mac.Ns1;

s0 = transpose(Mac.s(II));   R0 = R(II,:);
chi0 = Mac.chi;          Z0 = Z(II,:);
hs = min(s0(2:end)-s0(1:end-1))/2; hs = min(hs,2e-5);
hchi = min(chi0(2:end)-chi0(1:end-1))/2; hchi = min(hchi,1e-4);
s1 = s0 - hs;        s2 = s0 + hs;  
chi1 = chi0 - hchi;  chi2 = chi0 + hchi;  

% compute dR/ds using R(s,chi) and spline
R1 = spline(s0,R0',s1);  R2 = spline(s0,R0',s2);
dRds(II,:) = transpose((R2-R1)/hs/2); 

% compute dZ/ds using Z(s,chi) and spline
Z1 = spline(s0,Z0',s1);  Z2 = spline(s0,Z0',s2);
dZds(II,:) = transpose((Z2-Z1)/hs/2); 

% compute dR/dchi using R(s,chi) and spline
R1 = spline(chi0,R0,chi1);  R2 = spline(chi0,R0,chi2);
R1(:,1) = R1(:,end);  R2(:,end) = R2(:,1); 
dRdchi(II,:) = (R2-R1)/hchi/2; 

% compute dZ/dchi using Z(s,chi) and spline
Z1 = spline(chi0,Z0,chi1);  Z2 = spline(chi0,Z0,chi2);
Z1(:,1) = Z1(:,end);  Z2(:,end) = Z2(:,1); 
dZdchi(II,:) = (Z2-Z1)/hchi/2; 

% compute jacobian from (x,y,z) --> (s,chi,phi)
jacobian = (-dRdchi.*dZds + dRds.*dZdchi).*R;
jacobian(1,:) = jacobian(2,:);


% show jacobian
if Mac.plot_JACOB0
   figure(Mac.plot_JACOB0)
   %plot(Mac.chi,jacobian(Mac.Ns1,:),'r-'), hold on,
   N = round(Mac.Nchi/2);
   plot(R(:,1)*Mac.R0EXP,jacobian(:,1),'r-'), hold on,
   plot(R(:,N)*Mac.R0EXP,jacobian(:,N),'r-'), hold on,
   %plot(R(:,1)*Mac.R0EXP,dZds(:,1),'r-'), hold on,
   %plot(R(:,N)*Mac.R0EXP,dZds(:,N),'r-'), hold on,
   %plot(Mac.s,R(:,1)*Mac.R0EXP,'r-'), hold on,
   %plot(Mac.s,R(:,N)*Mac.R0EXP,'r-'), hold on,
   %x1=Mac.s; y1=R(:,1)*Mac.R0EXP;
   %x2=(x1(1:end-1)+x1(2:end))/2;
   %z2=diff(y1)./diff(x1);
   %plot(x2,z2,'b-+'), hold on,
end

% plot R,dRds,dRdchi
if 1==0
   [smin,II] = min(abs(Mac.s-1.206));
   figure(1)
   plot(Mac.chi,Z(II,:),'b--'), hold on,
   figure(2)
   plot(Mac.chi,dZds(II,:),'b--'), hold on,
   figure(3)
   plot(Mac.chi,dZdchi(II,:),'b--'), hold on,
end

if 1==0
   [smin,II] = min(abs(Mac.s-Mac.rw(1)));
   CHIS = Mac.CHIS*pi;
   Srs  = spline(Mac.chi,dRds(II,:),CHIS);
   Src  = spline(Mac.chi,dRdchi(II,:),CHIS);
   Szs  = spline(Mac.chi,dZds(II,:),CHIS);
   Szc  = spline(Mac.chi,dZdchi(II,:),CHIS);
   Sj   = spline(Mac.chi,jacobian(II,:),CHIS);
  
   for k=1:length(CHIS)
	    MAT_A   = [Srs(k) Src(k); Szs(k) Szc(k)]/Sj(k);
   end
end

if 1==1
%compute |grad s|_{m=0} inside plasma
G22   = dRdchi(1:Mac.Ns1,:).^2 + dZdchi(1:Mac.Ns1,:).^2;
G11U = G22.*(R(1:Mac.Ns1,:)./jacobian(1:Mac.Ns1,:)).^2;
grads = sum(sqrt(G11U),2)/(size(G11U,2)-1);
grads(1) = grads(2);
RES_gradS = [Mac.s(1:Mac.Ns1) grads];
save GRADSM0.OUT RES_gradS -ascii -double
end





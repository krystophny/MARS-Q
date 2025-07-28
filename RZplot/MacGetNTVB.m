% Convert the MARS-F computed B-field from
% EQARC coordinates (NEGP=-1, NER=1) to
% Hamada coordinates (NEGP=0, NER=0).
% The CHEASE computed coordinates with (NEGP=0, NER=0) 
% are not strictly Hamada coordinates (Jacobian=1), 
% with a difference in C(psi), which can be absorbed
% by re-defining s-coordinate (s will not be 1 at the plasma surface!).
% This should be taken into acount as soon as there is 
% partical differentiations along s, which is not the case for 
% the NTV calculations.
 
function BMH = MacGetNTVB(R,Z,RH,ZH,B)

global Mac

BMH = zeros(Mac.Ns1,Mac.Nm1);

%quick check
%R0Z0 = [R(1,:); Z(1,:); RH(1,:); ZH(1,:)]'

chi = Mac.chi(1:end-1);

for k=2:Mac.Ns1
    Zmin=min(Z(k,:));  Zmax=max(Z(k,:)); Z0=0.5*(Zmin+Zmax);
    Rmin=min(R(k,:));  Rmax=max(R(k,:)); R0=0.5*(Rmin+Rmax);
    T = atan2(Z(k,1:end-1)-Z0,R(k,1:end-1)-R0);
    [T,JJ] = sort(T);
    B0 = B(k,:);
    BB = B0(JJ);

    %transform B-field in Hamada coordinate chi-angle: Mac.chi
    TH = atan2(ZH(k,1:end-1)-Z0,RH(k,1:end-1)-R0);
    BH = spline(T,BB,TH);
    BH = [BH BH(1)];
    
    %get Gauss quadrature points for Fourier decomposition
    [z,w] = MacGaussQuad1D(4);
    x0 = (Mac.chi(1:end-1)+Mac.chi(2:end))*0.5;
    h2 = diff(Mac.chi)*0.5;
    xx = z'*h2 + ones(size(z'))*x0;
    wh = w'*h2;
    xx = xx(:);  wh=wh(:);

    %Fourier decompose BH in Hamada chi-angle
    mm = transpose(Mac.Mm(:));
    expmt = exp(-i*xx*mm)/(2*pi);
    BHH   = spline(Mac.chi,BH,xx').*wh';
    BHM   = BHH*expmt; 

    %store BHM to BMH
    BMH(k,:) = transpose(BHM(:));

    %plotting
    if k==Mac.Ns1 & Mac.plot_ChiHAM>0
      figure(Mac.plot_ChiHAM)
      subplot(2,1,1),
      plot(T,chi(JJ),'b-','LineWidth',2), hold on,
      [TH,JJ] = sort(TH);
      plot(TH,chi(JJ),'r-','LineWidth',2), hold on,
      xlabel('physical poloidal angle','FontSize',14)
      ylabel('equi-arc(blue) vs Hamada(red) poloidal angle','FontSize',14)

      subplot(2,1,2),
      plot(Mac.chi,BH,'r-','LineWidth',2), hold on,
      xlabel('Hamada poloidal angle','FontSize',14)
      ylabel('|B|','FontSize',14)
    end
end

BMH(1,:) = BMH(2,:);

if Mac.plot_BHAM > 0
  figure(Mac.plot_BHAM)
  subplot(3,1,1), plot(Mac.s(1:Mac.Ns1),real(BMH)), hold on,
                  ylabel('Re(B_{Hamada})','FontSize',14),
                  xlabel('s','FontSize',14),
  subplot(3,1,2), plot(Mac.s(1:Mac.Ns1),imag(BMH)), hold on,
                  ylabel('Im(B_{Hamada})','FontSize',14),
                  xlabel('s','FontSize',14),
  subplot(3,1,3), yy = abs(BMH).^2; 
                  plot(Mac.Mm,sqrt(sum(yy,1)/size(yy,1))), hold on,
                  ylabel('||B_{Hamada}||_{L_2}','FontSize',14),
                  xlabel('m','FontSize',14),
end

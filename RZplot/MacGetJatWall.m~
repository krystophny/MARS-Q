function [Rw,Zw,Jwr,Jwz,Jwphi] = MacGetJatWall(R,Z,Jr,Jz,Jphi)

global Mac

% find the index for wall surfaces
II = zeros(1,length(Mac.rw));
for j=1:length(Mac.rw)
    [smin,II(j)] = min(abs(Mac.s-Mac.rw(j)));
end

Tw = [];  Rw = [];  Zw = [];
Jwr = []; Jwz = []; Jwphi = [];
for j=1:length(Mac.rw)
    Zmin=min(Z(II(j),:));  Zmax=max(Z(II(j),:)); Z0=0.5*(Zmin+Zmax);
    Rmin=min(R(II(j),:));  Rmax=max(R(II(j),:)); R0=0.5*(Rmin+Rmax);
    T = atan2(Z(II(j),:)-Z0,R(II(j),:)-R0);
    I = find(T<0); T(I) = T(I) + 2*pi;  %for JET62653
    [T,JJ] = sort(T);
    Tw = [Tw T];
    Rw = [Rw R(II(j),JJ)];
    Zw = [Zw Z(II(j),JJ)];
    Jwr = [Jwr Jr(II(j),JJ)];
    Jwz = [Jwz Jz(II(j),JJ)];
    Jwphi = [Jwphi Jphi(II(j),JJ)];
end

data = [[Rw(:) Zw(:)]*Mac.R0EXP real(Jwr(:)) imag(Jwr(:)) real(Jwz(:)) imag(Jwz(:)) real(Jwphi(:)) imag(Jwphi(:))];

save MacDataJw data -ascii

if Mac.plot_JW3>0 
  figure(Mac.plot_JW3)
  N = Mac.Nchi;
  deltaw = 1.4103e-02; %[m], wall thickness 
  J00    = Mac.B0EXP/Mac.R0EXP/(4e-7*pi)*1e+6/deltaw;
  Jwphi  = Jwphi/J00;   %[MA/m]
  C      = max(abs(Jwphi));
  for j=1:length(Mac.rw)
    plot(Tw((j-1)*N+1:j*N)*180/pi,real(Jwphi((j-1)*N+1:j*N))/C,'b-','LineWidth',4), hold on,
    %plot(Tw((j-1)*N+1:j*N),imag(Jwr((j-1)*N+1:j*N)),'b-','LineWidth',2), hold on,
    %plot(Tw((j-1)*N+1:j*N),imag(Jwz((j-1)*N+1:j*N)),'k-','LineWidth',2), hold on,
    %plot(Tw((j-1)*N+1:j*N),imag(Jwphi((j-1)*N+1:j*N)),'r--','LineWidth',2), hold on,
  end
  xlabel('geometric poloidal angle [deg.]','FontSize',18,'FontWeight','Bold')
  ylabel('J_\phi^{surf}','FontSize',18,'FontWeight','Bold')
  ha = get(Mac.plot_JW3,'CurrentAxes');
  set(ha,'FontSize',18,'FontWeight','Bold')
  a=axis; axis([0 360 a(3) a(4)])

  %compute total wall current
  dl = sqrt((diff(Rw)).^2+(diff(Zw)).^2);
  Jwtot = sum((Jwphi(1:end-1)+Jwphi(2:end))/2.*dl)*Mac.R0EXP  %[MA]
  Lwtot = sum(dl)*Mac.R0EXP  %[m]
end

% plot Jw at the wall
if Mac.plot_JW>0
   figure(Mac.plot_JW)
   h = 3.0/Mac.Ns1;
   Jrr = real(Jwr); Jzr = real(Jwz);
   Jt = sqrt(Jrr.^2 + Jzr.^2);  
   %Jt = max(max(Jt));
   [II,JJ]=find(Jt==0.0);  Jt(II,JJ) = 1.0;
   R1 = Rw + h*Jrr./Jt;
   Z1 = Zw + h*Jzr./Jt;

   R2 = [Rw(:) R1(:)]';
   Z2 = [Zw(:) Z1(:)]';

   plot(Rw*Mac.R0EXP,Zw*Mac.R0EXP,'c-'), hold on,
   plot(R2*Mac.R0EXP,Z2*Mac.R0EXP,'b-'), hold on,
   axis equal
   xlabel('R [m]','FontSize',16)
   ylabel('Z [m]','FontSize',16)
end


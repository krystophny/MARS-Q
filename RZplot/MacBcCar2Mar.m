function  MacBcCar2Mar(filename)
%% Read Cariddi input coupling matrices for backward coupling scheme
%% and convert to pre-defined data format

format short e

global Fc 

eval(['load ' filename])

tauA = Fc.tauA;

for k=1:length(gamma2)
  if abs(gamma2(k))==0, g='0'; 
  %elseif abs(gamma2(k))==max(abs(gamma2)), g='inf'; 
  else g=num2str(gamma2(k)*tauA); end
  %end
  %g=num2str(gamma2(k)*tauA); 

  BU = [real(bn(:,:,k))' imag(bn(:,:,k))';
        real(bpol(:,:,k))' imag(bpol(:,:,k))';
        real(btor(:,:,k))' imag(btor(:,:,k))'];

  save tmp BU -ascii -double
  eval(['!mv tmp ' Fc.DIR 'BU' g '_ALL']);
end

Fc.gamma2 = gamma2;

kpsi = Fc.kpsi;

if kpsi > 0
   ii = [1:size(psi_x,1)];
   psi_x = psi_x(ii,:,:);
   psi_y = psi_y(ii,:,:);
   psi_z = psi_z(ii,:,:);
end

if kpsi == 1  %normal component
   Fc.psis = psi_x;
   N = size(psi_x,1);
   for k=1:N
       t = (k-1)*2*pi/N;
       Fc.psis(k,:,:) = psi_x(k,:,:)*cos(t) + psi_y(k,:,:)*sin(t);
   end
elseif kpsi == 2  %toroidal component
   Fc.psis = psi_x;
   N = size(psi_x,1);
   for k=1:N
       t = (k-1)*2*pi/N;
       Fc.psis(k,:,:) = psi_y(k,:,:)*cos(t) - psi_x(k,:,:)*sin(t);
   end
elseif kpsi == 3 %poloidal component
   Fc.psis = psi_z;
elseif kpsi == 4 %n=1 harmonic at theta=0, nomal component
   Fc.psis = sum(bn,1);
elseif kpsi == 5 %n=1 harmonic at theta=0, poloidal component
   Fc.psis = sum(bpol,1)*i;
elseif kpsi == 6 %n=1 harmonic at theta=0, toroidal component
   Fc.psis = sum(btor,1)*i;
end

% flux sensor for JET
if kpsi==7
   fac_sad  = 9.9144*1.018*16;
   fac_iron = 1.44;
   psi_flux = psi_flux*fac_sad*fac_iron;

   SX01 = psi_flux(1,:,:)*0.0524;
   SX14 = psi_flux(5,:,:)*(-0.0534);
   SY01 = psi_flux(2,:,:)*0.0536;
   SY14 = psi_flux(6,:,:)*(-0.0526);
   S101 = psi_flux(4,:,:)*0.0530;
   S114 = psi_flux(8,:,:)*(-0.0537);
   S501 = psi_flux(3,:,:)*0.0530;
   S514 = psi_flux(7,:,:)*(-0.0523);

	MHDF = (SX01 + SX14 - SY01 - SY14)/4;
	MHDG = (S101 + S114 - S501 - S514)/4;

   Fc.psis = zeros(1,size(psi_flux,2),size(psi_flux,3));
   Fc.psis(1,:,:) = MHDF-MHDG*i;
end

%compute the n=1 component of sensor signal using cubic spline 
if kpsi > 0 & kpsi < 4 
   n  = 1;
   N = size(Fc.psis,1);
   t0 = [0:N-1]*2*pi/N;
   M  = 201;
   h  = 2*pi/(M-1);
   t1 = [0:M-2]*h + h/2;
   g1 = exp(-i*n*t1);
   psis = Fc.psis;
   Fc.psis = sum(psis,1);
   for k=1:size(psis,2)
   for p=1:size(psis,3)
       y0 = psis(:,k,p);
       %t2 = [t0-2*pi t0 t0+2*pi];
       %y2 = [y0; y0; y0];
       t2 = [t0 2*pi];  
       y2 = [y0; y0(1)];
       y2r = real(y2);  y2i = imag(y2);
       a   = 1.5;
       y1 = interp1(t2,sign(y2r).*abs(y2r).^(1/a),t1,'linear') + ...
            i*interp1(t2,sign(y2i).*abs(y2i).^(1/a),t1,'linear');
       y1 = sign(real(y1)).*abs(real(y1)).^a + ...
            i*sign(imag(y1)).*abs(imag(y1)).^a;
       %y1 = spline(t2,y2,t1);
       Fc.psis(1,k,p) = sum(y1.*g1)*h/2/pi;

       %test plot
       if k==1 & p==1
          figure(11)
          plot(t0,real(y0),'ro',t1,real(y1),'b-')
          figure(12)
          plot(t0,imag(y0),'ro',t1,imag(y1),'b-')
       end
   end
   end
   %Q = size(psis,2);
   %p = exp(-i*n*2*pi/Q);
   %for k=2:Q
   %    Fc.psis(1,k,:) = Fc.psis(1,k-1,:)*p;
   %end
   if kpsi==2 | kpsi==3, Fc.psis = Fc.psis*i; end
end



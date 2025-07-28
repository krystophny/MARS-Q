% plot the plasma normal displacement

% input parameters
n    = 6;
phi0 =0/180*pi;  %displacement at the plane with phi=phi0
fac_scale = 100; %scaling factor

% read data
data = load('../MAST_ALL/ResLinResponse/XnSurf_27359_M59.asc');
R    = data(:,2);
Z    = data(:,3);
xir  = data(:,4);  
xii  = data(:,5);  

% compute normal displacement at the middle of each segment
YD   = diff(R) + diff(Z)*i;

xir = (xir(1:end-1)+xir(2:end))/2;
xii = (xii(1:end-1)+xii(2:end))/2;
Rx  = (R(1:end-1)+R(2:end))/2;
Zx  = (Z(1:end-1)+Z(2:end))/2;

Eta  = xir*cos(n*phi0) - xii*sin(n*phi0);
YN   = (Rx+Zx*i) + Eta.*YD./abs(YD)*(-i)*fac_scale;

% plotting...
figure(1)
plot(R,Z,'b--',real(YN),imag(YN),'r-','LineWidth',2), hold on,
axis equal
%axis([0.5 0.9 -1.35 -1.10])



% linear superposition of response data from different rows of RMP coils

SDIR1 = '/home/liuy/Work/IR_D/Database/DIII-D/157376/3D/n3/';
CP = 180;  %coil phasing [degree] between U and L coils

CC = exp(i*CP/180*pi);
SDIR2 = [SDIR1 'UL_' int2str(CP) '/'];
mkdir(SDIR2)

JU = load([SDIR1 'U/JPLASMA.OUT']);
JL = load([SDIR1 'L/JPLASMA.OUT']);

M   = round(JU(1,1));
II  = [0:3]*2;
JUU = (JU(M+2:end,II+1)+JU(M+2:end,II+2)*i);
JLL = (JL(M+2:end,II+1)+JL(M+2:end,II+2)*i);
JNN = JUU + JLL*CC;

JN  = [JU(1:M+1,:); real(JNN(:,1)) imag(JNN(:,1)) ...
                    real(JNN(:,2)) imag(JNN(:,2)) ...
                    real(JNN(:,3)) imag(JNN(:,3)) ...
                    real(JNN(:,4)) imag(JNN(:,4))];
save([SDIR2 'JPLASMA.OUT'],'JN','-ascii');

XU = load([SDIR1 'U/XPLASMA.OUT']);
XL = load([SDIR1 'L/XPLASMA.OUT']);

M   = round(XU(1,1));
N   = round(XU(1,2));
II  = [0:2]*2;
XUU = (XU(M+N+2:end,II+1)+XU(M+N+2:end,II+2)*i);
XLL = (XL(M+N+2:end,II+1)+XL(M+N+2:end,II+2)*i);
XNN = XUU + XLL*CC;

XN  = [XU(1:M+N+1,:); real(XNN(:,1)) imag(XNN(:,1)) ...
                      real(XNN(:,2)) imag(XNN(:,2)) ...  
                      real(XNN(:,3)) imag(XNN(:,3))];
save([SDIR2 'XPLASMA.OUT'],'XN','-ascii');








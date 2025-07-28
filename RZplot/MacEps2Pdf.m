%convert all eps-files from a folder to corresponding pdf-files

SDIRM = '/home/liuy/Codes/MarsQ/RZplot/';
%SDIRM = '/cscratch/liuy/Temp/';
%SDIRD = '/cscratch/liuy/Temp/';
%SDIRD = '/home/liuy/Work/DIII-D/D3D_QH_CETOP/Data/163520.2650/n1/';
SDIRD = '/home/liuy/Work/NSTX/NSTX132543/Work/';

cd(SDIRD);
FF = split(ls);
NF = length(FF)-1;

for kf=1:NF
    F = FF{kf};
    if strcmp(F(end-3:end),'.eps')
       eval(['!' SDIRM 'epstopdf ' F]);
    end   
end

cd(SDIRM)

 

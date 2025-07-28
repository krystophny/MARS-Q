close all
 
global SDIR rviw

kcase = 1;

if kcase == 1
%ITER
mm   = [-29:29];
rviw = 1.01;
SDIR0 = '/scratch/yliu/Iter/';

elseif kcase == 2
%CIRC
mm   = [-5:9];
rviw = 1.1;
SDIR0 = '/.automount/funsrv1/root/home/yliu/DataMarsf/';

elseif kcase == 3
%RFX
mm   = [-5:5];
rviw = 1.05;
SDIR0 = '/.automount/funsrv1/root/home/yliu/Rfx/';

elseif kcase == 4
%JT60-SA
mm   = [-9:29];
rviw = 1.05;
SDIR0 = '/scratch/yliu/';

elseif kcase == 5
%RFX tokamak configuration
mm   = [-3:9];
rviw = 1.01;
SDIR0 = '/.automount/funsrv1/root/home/yliu/Rfx/';
end

MGG = length(mm); 

MacFcCheckTot
MacFcCheckBiot
MacFcCheckRwm
MacFcCheckMat

eval(['!zip ITER_data' '_' g0 '.zip *.asc MacDataS'])
 
%eval(['save ITER_data' '_' g0 '.mat AtTot ApTot AtVac ApVac AnPls'])

%eval(['!rm -f *.asc MacDataS'])




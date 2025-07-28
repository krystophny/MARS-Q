close all, clear all,

global SDIR rviw

kcase = 1;

if kcase == 1
%ITER
mm   = [-29:29];
rviw = 1.01;
g0   = '1e-7';
SDIR0 = '/scratch/yliu/Iter/';

elseif kcase == 2
%CIRC
mm   = [-5:9];
rviw = 1.1;
g0   = '1e+6-2e-4i';
SDIR0 = '/.automount/funsrv1/root/home/yliu/DataMarsf/';

elseif kcase == 3
%RFX
mm   = [-5:5];
rviw = 1.05;
g0   = '1e-9';
SDIR0 = '/.automount/funsrv1/root/home/yliu/Rfx/';

elseif kcase == 4
%JT60-SA
mm   = [-9:29];
rviw = 1.05;
g0   = '1e-6';
SDIR0 = '/scratch/yliu/';

elseif kcase == 5
%RFX tokamak configuration
mm   = [-3:9];
rviw = 1.01;
g0   = '1e-7';
SDIR0 = '/.automount/funsrv1/root/home/yliu/Rfx/';

end

MGG = length(mm); 

MacFcGenMatSrfa
MacFcGenMatSbiot
MacFcGenMatSvac


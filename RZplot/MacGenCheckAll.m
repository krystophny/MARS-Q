close all
 
global SDIR rviw

kcase = 2;

if kcase == 1
%ITER
mm   = [-9:39];
rviw = 1.05;
SDIR = '/home/elfliu/FUSION/Scen4_V8/';

elseif kcase == 2
%CIRC
mm   = [-9:19];
rviw = 1.10;
SDIR = '/home/elfliu/FUSION/CariddiCirc/';
end

MGG = length(mm); 

MacGenCheckTot
MacGenCheckBiot
MacGenCheckRwm
MacGenCheckMatInv

!zip ITER_data.zip *.asc MacDataS
 





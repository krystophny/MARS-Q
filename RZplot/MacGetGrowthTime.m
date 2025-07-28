% find growth time of instability, for given: 
% - growth rate in MARS-F: gamma  
% - initial perturbation amplitude: perturb0
% - final perturbation amplitude: perturb1
% - Alfven time: TAUA

TAUA    = 2.9375057372504544E-007; %RE_ITER
PERTURB = [1.0220e-05; 1.0848e-05; 6.2150e-06]; %=1G

gamma    = [2.29434E-02; 1.00432E-02; 3.12114E-02]; %[68ms; 76ms; 86ms]
perturb0 = PERTURB*1.0;       %=1G
perturb1 = PERTURB*1.8854e+4; %=18.8kG
%perturb1 = [1; 1; 1];           %=B0EXP 

timeM = log(perturb1./perturb0)./gamma;
timeP = timeM*TAUA*1e+6;  %[microsecond]

growth_time = [timeM timeP]
 

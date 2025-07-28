MGG = 7; 
NRS = 1;
SDIR0 = '/home/yliu/DeltapAna/Work/';

B1PGG = zeros(MGG,MGG);  
B1AGG = zeros(MGG,NRS,MGG);
B1ATT = zeros(MGG,NRS);

for GenNo=1:MGG
  disp(['GenNo=',int2str(GenNo)])
  eval(['!cp ' SDIR0 'BPLASMA_Dp' int2str(GenNo) ' ' SDIR0 'BPLASMA.OUT']);
  MacMainDeltap_Dp
  B1PGG(:,GenNo) = BM1P;
  B1AGG(:,:,GenNo) = BM1A;
end

disp(['Condition number = ', num2str(cond(B1PGG))])  

C = inv(B1PGG);
b = zeros(MGG,1); 
b(4) = 1.0;
u = C*b;

% save BC for the full solution
res = [real(u) imag(u)];
save BNORM01 res -ascii -double
!mv BNORM01 ../DeltapAna/Work/BNORM01.IN

% save combined solution at all rational surfaces 
for GenNo=1:MGG
    B1ATT = B1ATT + B1AGG(:,:,GenNo)*u(GenNo);
end
save temp_B1ATT B1ATT

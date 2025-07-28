BGGn = zeros(MGG,MGG);  BGGt = BGGn;  BGGp = BGGn;

for GenNo=1:MGG
  disp(['GenNo=',int2str(GenNo)])
  SDIR = [SDIR0 'BPLASMA_Spn/'];
  MacMainITER_Srfa
  BGGn(:,GenNo) = Bnm(:);
  MacMainITER_Sbiot
  BGGnB(:,GenNo) = Bnm(:);
end

disp(['Condition number = ', num2str(cond(BGGn))])  %=1.2561e+01

C = inv(BGGn);
At = BGGnB*C;  AnPls = At; At = [real(At) imag(At)];

save ResMatrixBiot_n.asc At -ascii -double
 

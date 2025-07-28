BGGn = zeros(MGG,MGG);  BGGt = BGGn;  BGGp = BGGn;

for GenNo=1:MGG
  disp(['GenNo=',int2str(GenNo)])
  SDIR = [SDIR0 'BPLASMA_Spn/'];
  MacMainITER_Srfa
  BGGn(:,GenNo) = Bnm(:);
  BGGt(:,GenNo) = Btm(:);
  BGGp(:,GenNo) = Bpm(:);
end

disp(['Condition number = ', num2str(cond(BGGn))])  %=1.2561e+01

C = inv(BGGn);
At = BGGt*C;  AtTot=At; At = [real(At) imag(At)];
Ap = BGGp*C;  ApTot=Ap; Ap = [real(Ap) imag(Ap)];

save ResMatrixTot_t.asc At -ascii -double
save ResMatrixTot_p.asc Ap -ascii -double


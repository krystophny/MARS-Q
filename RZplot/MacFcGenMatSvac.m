BGGn = zeros(MGG,MGG);  BGGt = BGGn;  BGGp = BGGn;

for GenNo=1:MGG
  disp(['GenNo=',int2str(GenNo)])
  SDIR = [SDIR0 'BPLASMA_Svn/'];
  MacMainITER_Svac
  BGGn(:,GenNo) = Bnm(:);
  BGGt(:,GenNo) = Btm(:);
  BGGp(:,GenNo) = Bpm(:);
end

disp(['Condition number = ', num2str(cond(BGGn))])  %=1.2561e+01

C = inv(BGGn);
At = BGGt*C;  AtVac = At;  At = [real(At) imag(At)];
Ap = BGGp*C;  ApVac = Ap; Ap = [real(Ap) imag(Ap)];

save ResMatrixVac_t.asc At -ascii -double
save ResMatrixVac_p.asc Ap -ascii -double
 

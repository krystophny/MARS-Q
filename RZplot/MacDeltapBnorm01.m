% generating BNORM01 for BC in Deltap computations

mm = [-5:10]';
M  = length(mm);

for k = 1:M
  Bnm = zeros(size(mm));
  Bnm(k) = 1.0;
  Bn = [real(Bnm) imag(Bnm)];

  scom = ['save BNORM01d' num2str(k) ' Bn -ascii'];
  eval(scom)
end

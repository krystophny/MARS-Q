%generate BNORM01 files for strong coupling scheme

M1 = -5; 
M2 = 9;

mm = M1:M2;
Bn0 = zeros(length(mm),2);

for k=1:length(mm)
  m = k - 1 + M1;
  Bn = Bn0;
  Bn(k,1) = 1.0;
  
  scom = ['save BNORM01o' num2str(k) ' Bn -ascii'];
  eval(scom)
end
  

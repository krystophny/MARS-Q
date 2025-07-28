global chis Bss1 Bs10 Bs20 knorm chis0

chis0 = 0.0;

!cp /home/elfliu/FUSION/JetSymData/PLASMAdata08/BPLASMA_SCW /home/elfliu/FUSION/JetSymData/PLASMAdata08/BPLASMA_TMP
knorm = 1;
MacMainJetSym; Bv = Bss1; Bv1 = Bs10;
knorm = 0;

MM = 23;

B = zeros(MM,length(Bv));  B1 = zeros(MM,1);  B2 = B1;

!cp /home/elfliu/FUSION/JetSymData/PLASMAdata03/BPLASMA_RWM /home/elfliu/FUSION/JetSymData/PLASMAdata08/BPLASMA_TMP
MacMainJetSym; B(1,:) = Bss1; B1(1) = Bs10; B2(1) = Bs20;

ss1 = ['02'; '05'; '06'; '09'; '11'; '12'; '18'; '21'; '22'];
ss2 = ['01'; '03'; '04'; '07'; '08'; '10'; '13'; '14'; '15'; '16'; '17'; '19'; '20';];
ss3 = ['01'; '02'; '03'; '04'; '05'; '06'; '07'; '08'; '09'; '10'; '11'; '12'; '13'; '14'; '15'; '16'; '17'; '18'; '19'; '20'; '21'; '22'];
ss =ss3;

for k=1:MM-1
  eval(['!cp /home/elfliu/FUSION/JetSymData/PLASMAdata03/BPLASMA_' ss(k,:) ' /home/elfliu/FUSION/JetSymData/PLASMAdata08/BPLASMA_TMP']);
  MacMainJetSym; 
  B(k+1,:) = Bss1; B1(k+1) = Bs10; B2(k+1) = Bs20;
end

A = zeros(MM,MM); b = zeros(MM,1);
chih = chis(2)-chis(1);  ff = exp(-chis.^2/pi^2);
for k=1:MM
  b(k) = (sum(ff.*Bv.*conj(B(k,:)))-ff(1)*Bv(1)*conj(B(k,1))*0.5-ff(end)*Bv(end)*conj(B(k,end))*0.5)*chih;
  for m=1:k
    A(k,m) = (sum(ff.*B(m,:).*conj(B(k,:)))-ff(1)*B(m,1)*conj(B(k,1))*0.5-ff(end)*B(m,end)*conj(B(k,end))*0.5)*chih;
  end
end

for k=1:MM
  for m=1:k-1
    A(m,k) = conj(A(k,m));  
  end
end

c = A\b;


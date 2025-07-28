global chis Bss1 Bs10 Bs20 Bs30 knorm chis0

knorm = 0; chis0 = 0.2;
MM = 50;
Br = zeros(MM,1);  Bpi = Br; Bpe = Br;

for k=1:MM
  eval(['!cp /home/elf/elfliu/FUSION/JET/SHOT40542/Results_2/BPLASMAF_' num2str(k) ' /home/elf/elfliu/FUSION/JET/SHOT40542/Results_2/BPLASMA_TMP']);
  MacMainJetSym; 
  Br(k) = Bs10; Bpi(k) = Bs20*1i; Bpe(k) = Bs30*1i;
end

Ball = [real(Br) imag(Br) real(Bpe) imag(Bpe) real(Bpi) imag(Bpi)];
 
for k=1:MM
  eval(['!cp /home/elf/elfliu/FUSION/JET/SHOT40542/Results_2/BPLASMAG_' num2str(k) ' /home/elf/elfliu/FUSION/JET/SHOT40542/Results_2/BPLASMA_TMP']);
  MacMainJetSym; 
  Br(k) = Bs10; Bpi(k) = Bs20*1i; Bpe(k) = Bs30*1i;
end

Ball = [Ball; real(Br) imag(Br) real(Bpe) imag(Bpe) real(Bpi) imag(Bpi)];
 
save /home/elf/elfliu/FUSION/JET/SHOT40542/Results_2/FullModelData Ball -ascii


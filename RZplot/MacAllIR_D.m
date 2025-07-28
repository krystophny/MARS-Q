global SDIR

kkk = 4;  %1: create new folders
          %2: after vacuum runs
          %3: after plasma response runs
          %4: save PROFEQ.OUT and RMZM_F.OUT
nnn = 4;

SS1 = 'ABCD';
SS2 = 'ULM';

for k=1:length(SS1)
for l=1:length(SS2)  
    if kkk==1
       mkdir(['~/Work/IR_D/' SS1(k) '/n' int2str(nnn) '/' SS2(l)])
    end

    if kkk==2 
       if k==1 
          SDIR = ['/cscratch/liuy/WorkIR_D/' SS1(k) SS2(l) '/'];
          MacMainIR_D
       end
       copyfile('dBnormal.txt',['~/Work/IR_D/' SS1(k) '/n' int2str(nnn) '/' SS2(l) '/dBnormal.txt'],'f')
     end

     if kkk==3
        copyfile(['/cscratch/liuy/WorkIR_D/' SS1(k) SS2(l) '/BPLASMA.OUT'],['~/Work/IR_D/' SS1(k) '/n' int2str(nnn) '/' SS2(l) '/BPLASMA.OUT'],'f')
        copyfile(['/cscratch/liuy/WorkIR_D/' SS1(k) SS2(l) '/JPLASMA.OUT'],['~/Work/IR_D/' SS1(k) '/n' int2str(nnn) '/' SS2(l) '/JPLASMA.OUT'],'f')
        copyfile(['/cscratch/liuy/WorkIR_D/' SS1(k) SS2(l) '/XPLASMA.OUT'],['~/Work/IR_D/' SS1(k) '/n' int2str(nnn) '/' SS2(l) '/XPLASMA.OUT'],'f')
        copyfile(['/cscratch/liuy/WorkIR_D/' SS1(k) SS2(l) '/PPLASMA.OUT'],['~/Work/IR_D/' SS1(k) '/n' int2str(nnn) '/' SS2(l) '/PPLASMA.OUT'],'f')
     end
end
end

if kkk==4
   for k=1:length(SS1)
       copyfile(['/cscratch/liuy/WorkIR_D/' SS1(k) 'U/PROFEQ.OUT'],['~/Work/IR_D/' SS1(k) '/n' int2str(nnn) '/PROFEQ.OUT'],'f')
       copyfile(['/cscratch/liuy/WorkIR_D/' SS1(k) 'U/RMZM_F.OUT'],['~/Work/IR_D/' SS1(k) '/n' int2str(nnn) '/RMZM_F.OUT'],'f')
   end
end


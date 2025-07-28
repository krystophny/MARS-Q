MTs   = '0abcdefgh';
%MTdir = '/.automount/funsrv1/root/home/yliu/MAST020333/New/';
MTdir = '/.automount/funsrv1/root/home/yliu/Scen2_V01/';
MTres = [];

for k=1:length(MTs);
    MTc = MTs(k);
    eval(['!cp /scratch/yliu/Iter/n=4/BPLASMA' MTc ' ' MTdir 'BPLASMA']);
    MacRfaCtBn2
    MTres = [MTres BnPEST_RS(3,:)'];
    %eval(['!cp ' MTdir 'BnMat.mat ' MTdir 'BnMat_Q' MTc]);
end

MTres = [mq' MTres];
save DataB1mrot MTres -ascii
eval(['!mv DataB1mrot ' MTdir 'DataB1mrot_n4']);

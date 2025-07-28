% produce data for 2D plot of Jphi

MTs   = 'abcdefgh';
%MTdir = '/.automount/funsrv1/root/home/yliu/MAST020333/New/';
MTdir = '/.automount/funsrv1/root/home/yliu/Scen2_V01/';
MTj = [];

for k=1:length(MTs);
    MTc = MTs(k);
    %eval(['!cp /scratch/yliu/MAST/JPLASMA_P' MTc ' ' MTdir 'JPLASMA']);
    eval(['!cp /scratch/yliu/Iter/n=4/JPLASMA' MTc ' ' MTdir 'JPLASMA']);
    MacMainITER_15MA
    [X,II]= min(abs(Mac.chi));
    MTj   = [MTj Jphi(1:Mac.Ns1,II)];
end

MTj = [Mac.s(1:Mac.Ns1) MTj];
MTj = [MTj(:,1) real(MTj(:,2:end)) imag(MTj(:,2:end))];

save DataJphisrot MTj -ascii
eval(['!mv DataJphisrot ' MTdir 'DataJphisrot_n4']);

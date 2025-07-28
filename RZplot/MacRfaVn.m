% produce data for 2D plot of Vn

MTs   = 'abcdefgh';
%MTdir = '/.automount/funsrv1/root/home/yliu/MAST020333/New/';
MTdir = '/.automount/funsrv1/root/home/yliu/Scen2_V01/';
MTv = [];

I0 = 7;

for k=1:length(MTs);
    MTc = MTs(k);
    eval(['!cp /scratch/yliu/Iter/n=4/VPLASMA' MTc ' ' MTdir 'VPLASMA']);
    %MacMainMAST020333
    MacMainITER_15MA
    MTv   = [MTv transpose(Vn(Mac.Ns1-I0,:))];
end

Rs = R(Mac.Ns1-I0,:);
Zs = Z(Mac.Ns1-I0,:);
Rc = (min(Rs)+max(Rs))/2;
Zc = (min(Zs)+max(Zs))/2;
Tg = atan2(Zs-Zc,Rs-Rc);
MTv = [Tg' MTv];

MTv = [MTv(:,1) real(MTv(:,2:end)) imag(MTv(:,2:end))];

save DataVntrot MTv -ascii
eval(['!mv DataVntrot ' MTdir 'DataVntrot_n4']);

% setup runs for coil phasing scan
% for NTV computations in ITER

kcase = 3;

dataFEEDI_n1 = [
8.1860e-05   5.8507e-05   8.4844e-05
5.5604e-05   3.9741e-05   5.7631e-05
2.7802e-05   1.9871e-05   2.8815e-05
	     ];
dataFEEDI_n2 = [
7.9079e-05   5.7539e-05   8.1860e-05
5.3715e-05   3.9083e-05   5.5604e-05
2.6857e-05   1.9542e-05   2.7802e-05
		];
dataFEEDI_n3 = [
7.4823e-05   5.5940e-05   7.7018e-05
5.0824e-05   3.7997e-05   5.2315e-05
2.5412e-05   1.8999e-05   2.6158e-05
		];
dataFEEDI_n4 = [
6.8979e-05   5.3755e-05   7.0533e-05
4.6854e-05   3.6513e-05   4.7910e-05
2.3427e-05   1.8257e-05   2.3955e-05
		];

dataFEEDI = dataFEEDI_n3;

if kcase==1
   SDIR = '/cscratch/liuy/WorkIR_D/AT/';
end
if kcase==2
   SDIR = '/cscratch/liuy/WorkIR_D/BT/';
end
if kcase==3
   SDIR = '/cscratch/liuy/WorkIR_D/CT/';
end
   
FI = dataFEEDI(kcase,:);

for k=1:8
for l=1:8
    FIN = [FI(1)*exp(i*(k-1)*2*pi/8) FI(2) FI(3)*exp(i*(l-1)*2*pi/8)];
    sf  = ['RmarsQ_' int2str(k) int2str(l)]; 
    copyfile([SDIR 'RmarsQ'],[SDIR sf],'f');
    s1 = sprintf('%11.4e',real(FIN(1)));
    s2 = sprintf('%11.4e',imag(FIN(1)));
    s3 = sprintf('%11.4e',real(FIN(2)));
    s4 = sprintf('%11.4e',imag(FIN(2)));
    s5 = sprintf('%11.4e',real(FIN(3)));
    s6 = sprintf('%11.4e',imag(FIN(3)));
    ss = ['FEEDI  = (' s1 ',' s2 '),(' s3 ',' s4 '),(' s5 ',' s6 '),'];
    eval(['!sed -i s/''FEEDI  = (1.0,0.0),''/''' ss '''/g ' SDIR sf]);

    if 1==0
    fid = fopen('tmp','w');
    fprintf(fid,'%s\n',['./' sf]);
    fprintf(fid,'%s\n',['grep TORQUE log_mars > TORQUE_' int2str(k) int2str(l)]);
    fclose(fid); 
    eval(['!cat tmp >> ' SDIR 'all.bat']);  
    end 
end
end

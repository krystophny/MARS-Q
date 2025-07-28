% generate input radial profiles for MARS-*
% from EFIT p-file

R0EXP = 1.6995;
%B0EXP = 1.9169;  %184003.2750
B0EXP = 1.9140;  %184003.3880
%B0EXP = 1.9098;  %184003.2740
%psipp   = [0.9361579      0.9469605     0.9675397     0.9765310      0.9912874]; %2740
psipp   = [0.9370174      0.9455787     0.9686599     0.9767946      1]; %3880

SDIR   = '/home/liuy/Work/D3D184003/';
shot   = '184003.3880';
kprint = 0;

% read data from p-file
pname = [SDIR 'Expt/p' shot];

fid = fopen(pname,'r');
N   = fscanf(fid,'%i');
fclose(fid);
fid = fopen(pname,'r');
psi = zeros(N,1);

ne  = psi;  fgetl(fid);  %[10^20/m^3]
for k=1:N
    a = str2num(fgetl(fid));
    psi(k) = a(1);
    ne(k)  = a(2);
end

te  = psi;  fgetl(fid);  %[KeV]
for k=1:N
    a = str2num(fgetl(fid));
    te(k)  = a(2);
end

ni  = psi;  fgetl(fid);  %[10^20/m^3]
for k=1:N
    a = str2num(fgetl(fid));
    ni(k)  = a(2);
end

ti  = psi;  fgetl(fid);  %[KeV]
for k=1:N
    a = str2num(fgetl(fid));
    ti(k)  = a(2);
end

nb  = psi;  fgetl(fid);  %[10^20/m^3]
for k=1:N
    a = str2num(fgetl(fid));
    nb(k)  = a(2);
end

pb  = psi;  fgetl(fid);  %[Pa]
for k=1:N
    a = str2num(fgetl(fid));
    pb(k)  = a(2);
end

ptot  = psi;  fgetl(fid);  %[KPa]
for k=1:N
    a = str2num(fgetl(fid));
    ptot(k)  = a(2);
end

omeg  = psi;  fgetl(fid);  %[kRad/s]
for k=1:N
    a = str2num(fgetl(fid));
    omeg(k)  = a(2);
end

omegp  = psi;  fgetl(fid);  %[kRad/s]
for k=1:N
    a = str2num(fgetl(fid));
    omegp(k)  = a(2);
end

omgvb  = psi;  fgetl(fid);  %[kRad/s]
for k=1:N
    a = str2num(fgetl(fid));
    omgvb(k)  = a(2);
end

omgpp  = psi;  fgetl(fid);  %[kRad/s]
for k=1:N
    a = str2num(fgetl(fid));
    omgpp(k)  = a(2);
end

omgeb  = psi;  fgetl(fid);  %[kRad/s]
for k=1:N
    a = str2num(fgetl(fid));
    omgeb(k)  = a(2);
end

er  = psi;  fgetl(fid);  %[kV/m]
for k=1:N
    a = str2num(fgetl(fid));
    er(k)  = a(2);
end

fclose(fid);

s   = sqrt(psi);

hf=figure(1);
plot(psi,ne,'b-',psi,ni,'r--',psi,nb,'k:','LineWidth',3), hold on,
xlabel('\psi_n','FontSize',18,'FontWeight','Bold')
ylabel('number density [10^{20}/m^3]','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('N_e','N_i','N_{EP}')
if kprint>0, print(1,[SDIR 'Resu/ProfDen_D3D' shot '.eps'],'-depsc'), end

hf=figure(2);
plot(psi,te,'b-',psi,ti,'r--','LineWidth',3), hold on,
xlabel('\psi_n','FontSize',18,'FontWeight','Bold')
ylabel('temperature [keV]','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('T_e','T_i')
if kprint>0, print(2,[SDIR 'Resu/ProfT_D3D' shot '.eps'],'-depsc'), end

hf=figure(3);
omgtt = omgeb - omgpp;
plot(psi,omeg,'b-',psi,omgeb,'r:','LineWidth',3), hold on,
%plot(psi,omgpp,'k:',psi,omgtt,'g-',psi,omgpp,'k:','LineWidth',3), hold on,
a=axis;
plot([0 1],[0 0],'k-')
for k=1:length(psipp)
    plot([psipp(k) psipp(k)],[a(3) a(4)],'k--')
end
xlabel('\psi_n','FontSize',18,'FontWeight','Bold')
ylabel('frequencies [krad/s]','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('\Omega','\Omega_{ExB}')
%legend('\Omega','\Omega_{ExB}','\Omega_*','\Omega_{ExB}-\Omega_*')
if kprint>0, print(3,[SDIR 'Resu/ProfFreq_D3D' shot '.eps'],'-depsc'), end

hf=figure(4);
plot(psi,ptot,'b-',psi,pb/1000,'k:','LineWidth',3), hold on,
xlabel('\psi_n','FontSize',18,'FontWeight','Bold')
ylabel('pressure [kPa]','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('P_{tot}','P_{EP}')
if kprint>0, print(4,[SDIR 'Resu/ProfP_D3D' shot '.eps'],'-depsc'), end

format short e
res = [ne(1)*1e+20 ti(1)*1e+3 te(1)*1e+3 omeg(1)*1e+3 omgeb(1)*1e+3]

% saving to PROF*.IN input data for MARS-Q
fid = fopen([SDIR 'Data/PROFDEN_D3D' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; ne'/ne(1)]);
fclose(fid);

fid = fopen([SDIR 'Data/PROFROT_D3D' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; omeg'*1000]);
fclose(fid);

fid = fopen([SDIR 'Data/PROFWE_D3D' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; omgeb'*1000]);
fclose(fid);

fid = fopen([SDIR 'Data/PROFTI_D3D' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; ti'*1000]);
fclose(fid);

fid = fopen([SDIR 'Data/PROFTE_D3D' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; te'*1000]);
fclose(fid);

fid = fopen([SDIR 'Data/PROFDA_D3D' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; nb'./ne']);
fclose(fid);

fid = fopen([SDIR 'Data/PROFPA_D3D' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; pb'./ptot'/1000]);
fclose(fid);

 


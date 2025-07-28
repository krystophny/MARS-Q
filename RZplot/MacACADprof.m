function RES=MacACADprof(SDIR1,SDIR2,pname)

% generate input radial profiles for MARS-Q
% from EFIT p-file

fid = fopen([SDIR1 pname],'r');
N   = fscanf(fid,'%i');
fclose(fid);

fid = fopen([SDIR1 pname],'r');
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

pb  = psi;  fgetl(fid);  %[KPa]
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

if 1==0
hf=figure(1);
plot(psi,ne,'b-',psi,ni,'r--','LineWidth',2), hold on,
xlabel('\psi_n','FontSize',16,'FontWeight','Bold')
ylabel('Ne,Ni [10^{20}/m^3]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('Ne','Ni')

hf=figure(2);
plot(psi,te,'b-',psi,ti,'r--','LineWidth',2), hold on,
xlabel('\psi_n','FontSize',16,'FontWeight','Bold')
ylabel('Te,Ti [keV]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('Te','Ti')

hf=figure(3);
plot(psi,omeg,'b-',psi,omgeb,'r--','LineWidth',2), hold on,
xlabel('\psi_n','FontSize',16,'FontWeight','Bold')
ylabel('\Omega [krad/s]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
legend('\Omega','\Omega_{ExB}')

hf=figure(4);
plot(psi,ptot,'b-','LineWidth',2), hold on,
xlabel('\psi_n','FontSize',16,'FontWeight','Bold')
ylabel('Ptot [kPa]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
end

RES  = [ne(1)*1e+20 ti(1)*1e+3 te(1)*1e+3 omeg(1)*1e+3 omgeb(1)*1e+3];
%shot = pname(2:end);
shot = 'SAVE';

% saving to PROF*.IN input data for MARS-Q
fid = fopen([SDIR2 'PROFDEN_' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; ne'/ne(1)]);
fclose(fid);

fid = fopen([SDIR2 'PROFROT_' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; omeg'*1000]);
fclose(fid);

fid = fopen([SDIR2 'PROFWE_' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; omgeb'*1000]);
fclose(fid);

fid = fopen([SDIR2 'PROFTI_' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; ti'*1000]);
fclose(fid);

fid = fopen([SDIR2 'PROFTE_' shot '.IN'],'w');
fprintf(fid,'%4i  %4i\n',[length(s); 1]);
fprintf(fid,'%13.5e  %13.5e\n',[s'; te'*1000]);
fclose(fid);

 


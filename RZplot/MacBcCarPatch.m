function  MacBcCarPatch(filename)
%% Read Cariddi input coupling matrices for backward coupling scheme
%% for single coils response
%% and generate an array of coils.

eval(['load ' filename])

N = size(bn); M = N(1); K = N(2); J = N(3);
N = size(psi_x); L = N(1);

if K==1

Kn = L;

bnn = zeros(M,Kn,J); bpoln = bnn; btorn = bnn;
psi_xn = zeros(L,Kn,J); psi_yn = psi_xn; psi_xn = psi_xn;

psi_x2 = [psi_x; psi_x];
psi_y2 = [psi_y; psi_y];
psi_z2 = [psi_z; psi_z];

for k=1:Kn
    p = exp(i*2*pi*(k-1)/Kn);
    bnn(:,k,:) = bn(:,1,:)*p;
    bpoln(:,k,:) = bpol(:,1,:)*p;
    btorn(:,k,:) = btor(:,1,:)*p;
    II = [1:L] + k - 1;
    psi_xn(:,k,:) = psi_x2(II,1,:);
    psi_yn(:,k,:) = psi_y2(II,1,:);
    psi_zn(:,k,:) = psi_z2(II,1,:);
end;

bn = bnn;  bpol = bpoln; btor = btorn;
psi_x = psi_xn; psi_y = psi_yn; psi_z = psi_zn;

save tmp bn bpol btor psi_x psi_y psi_z gamma2
eval(['!mv tmp.mat ' filename])
 
end    

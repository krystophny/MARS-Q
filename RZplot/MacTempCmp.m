% compare BPLASMA_MARSF data (old vs. new)
% for Heinke

M1 = 85;
M2 = 85;
N1 = 786;
N2 = 837;

SDIR1 = '/cscratch/liuy/WorkTEMP/Old/';
SDIR2 = '/cscratch/liuy/WorkTEMP/New/';

s1 = load([SDIR1 'S']);
s2 = load([SDIR2 'S']);
b1 = load([SDIR1 'B']);
b2 = load([SDIR2 'B']);

s1 = s1(:,1);
s2 = s2(:,1);

b11 = reshape(b1(:,1),N1,M1);
b21 = reshape(b2(:,1),N2,M2);

figure(1)
mm = 4:4;
plot(s1,b11(:,mm+42+1),'r-','LineWidth',2), hold on
plot(s2,b21(:,mm+42+1),'b-','LineWidth',2)
a=axis; axis([0 1 a(3) a(4)])

mesh1 = load([SDIR1 'SCHIMESH_RECTRZ.IN']);
mesh2 = load([SDIR2 'SCHIMESH_RECTRZ.IN']);

NR = mesh1(1,1); NZ = mesh1(1,2);
z1 = reshape(mesh1(2:end,2),NZ,NR);
c1 = reshape(mesh1(2:end,4),NZ,NR);

NR = mesh2(1,1); NZ = mesh2(1,2);
z2 = reshape(mesh2(2:end,2),NZ,NR);
c2 = reshape(mesh2(2:end,4),NZ,NR);

figure(2)
II = 100;
plot(z1(:,II),c1(:,II),'r-','LineWidth',2), hold on
plot(z2(:,II),c2(:,II),'b-','LineWidth',2), hold on


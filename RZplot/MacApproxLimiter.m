% plot boundary and limiter shapes 
% approximate rectangular domain boundary of g-file by a smooth shape
% to be used in EXPEQ as a new wall surface
% this allows better enclusion of the limiter surface by the computational domain in vacuum
% from LIMITER.OUT produced by MARS-F (REORBIT module)

function MacApproxLimiter(kfig,R0EXP)

%SDIR = '/home/liuy/Work/MAST-U45163/Data/';
SDIR = '/home/liuy/Work/KSTAR/Data/30306.7850_gpec/';

d = load([SDIR 'LIMITER.OUT']);

RLFT=d(1,1); RLEN=d(1,2); RRGT=RLFT+RLEN;
ZLFT=d(2,1); ZLEN=d(2,2); ZRGT=ZLFT+ZLEN;

hf = figure(kfig);
x = [RLFT RRGT RRGT RLFT RLFT];
y = [ZLFT ZLFT ZRGT ZRGT ZLFT];

N = 76;
a1 = ones(1,N);             a1 = a1(2:N-1);
x1 = linspace(RLFT,RRGT,N); x1 = x1(2:N-1);
y1 = linspace(ZLFT,ZRGT,N); y1 = y1(2:N-1);

xx = [RLFT x1 RRGT RRGT*a1 RRGT fliplr(x1) RLFT RLFT*a1 RLFT];
yy = [ZLFT ZLFT*a1 ZLFT y1 ZRGT ZRGT*a1 ZRGT fliplr(y1) ZLFT]; 


xa = max(x)-min(x);
ya = max(y)-min(y);

x2 = xx/xa;
y2 = yy/ya;

x0 = (max(x2)+min(x2))/2;
y0 = (max(y2)+min(y2))/2;

rr = ((x2-x0).^2+(y2-y0).^2).^0.5;
rmin = min(rr);
rmax = max(rr);

x3 = x0 + (x2-x0)*rmin./rr;
y3 = y0 + (y2-y0)*rmin./rr;

x4 = x3*xa;
y4 = y3*ya;

plot(x4,y4,'r-+')
axis equal

%R0EXP = 0.8;
res = [x4(:) y4(:)]/R0EXP;
save BOUNDW res -ascii



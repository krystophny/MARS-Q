d = load('/cscratch/liuy/ZhaoYF/TEMP03');

t = d(:,1);

I = find(t>300);
t = t(I);
d = d(I,:);
s = d(:,2);

for k=2:size(d,2)
    figure(k)
    plot(t,d(:,k),'b-+')
end


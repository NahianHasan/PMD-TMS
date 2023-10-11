clear all
A=zeros([1077 14994]);
B=zeros([14994 1077]);
for i=1:1077
load(strcat('./fullmatrixinfo/soln2col',num2str(i),'.mat'));
A(:,i)=Efield;
end
for i=1:1077
load(strcat('./fullmatrixinfo/soln2row',num2str(i),'.mat'));
B(:,i)=Hobs;
end
[u,v,erri,capI,capJ]=ACA(B,200);
nori=norm(B);
for i=1:200;
X=u(:,1:i)*v(:,1:i)';
trueerr(i)=norm(B-X)/nori;
end
%%
semilogy(erri)

hold on
semilogy(trueerr,'r')
s=svd(B);

semilogy((2:1077)'.*s(2:end)/s(1),'green')

semilogy(s(2:end)/s(1),'black')
function [u,v,erri,capI,capJ]=ACA(A,Maxmodes);
%%%%%%%%begin ACA
capJ=ones([Maxmodes 1]);
norhist=zeros([Maxmodes 1]);
capI=ones([Maxmodes+1 1]);
[ncols nrows]=size(A);
v=zeros([nrows Maxmodes]);
u=zeros([ncols Maxmodes]);
erri=zeros([Maxmodes 1]);
for i=1:Maxmodes
%compute row
v(:,i)=A(capI(i),:)';
if i~=1
    v(:,i)=v(:,i)-(u(capI(i),1:i-1)*v(:,1:i-1)')';
end
%define J
 C = setdiff(1:nrows,capJ);
pp=find(abs(v(C,i))==max(abs(v(C,i))));
capJ(i)=C(pp(1));
v(:,i)=v(:,i)/v(capJ(i),i);
%compute column
u(:,i)=A(:,capJ(i));
if i~=1
u(:,i)=u(:,i)-u(:,1:i-1)*v(capJ(i),1:i-1)';
end
ukvk=norm(u(:,i))*norm(v(:,i));
norhist(i)=ukvk^2;
if i~=1
norhist(i)=norhist(i)+norhist(i-1)+2*sum(abs(u(:,1:i-1)'*u(:,i)).*abs(v(:,1:i-1)'*v(:,i)));
end
    %compute next index
    
 C = setdiff(1:ncols,capI);
pp=find(abs(u(C,i))==max(abs(u(C,i))));
capI(i+1)=C(pp(1));
erri(i)=ukvk/sqrt(norhist(i));

end
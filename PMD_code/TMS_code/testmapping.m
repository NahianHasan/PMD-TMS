load example1.mat
rs=rs';ks=ks';
that1(3,1:181)=that(3,1);
that1(2,1:181)=that(1,1);
that1(1,1:181)=that(2,1);
clear robs
[X,Y,Z]=ndgrid(-.005:.0003:.005,.06-.005:.0003:.06,-.005:.0003:.005);
robs(1,:)=X(:)';
robs(2,:)=Y(:)';
robs(3,:)=Z(:)';
[te2te]=gente2te(te2p);
[teid,bari]=pointlocation_c(robs,te2p-1,p,te2te,4);
teid=unique(teid(:));
that1=zeros([3,numel(teid)]);
that1(3,:)=1;

[rv,jv]=runcoderecip(te2p,p,conductivity,teid(:)-1,that1,3);
%%
N=[7 7 2];
Nj=numel(rs)/3;
rs2=zeros([3 Nj 360]);
ks2=zeros([3 Nj 360]);

for i=1:360
    phi=i*pi/180;
rs2(3,:,i)=rs(3,:)-.09;
ks2(:,:,i)=ks;
rs2(1,:,i)=rs(1,:)*cos(phi)+rs(2,:)*sin(phi);
rs2(2,:,i)=-rs(1,:)*sin(phi)+rs(2,:)*cos(phi);
end
% tic
% [Hobs]=computeHprimary(rv13,jv13,numel(rv13)/3,reshape(rs2,[3,360*Nj]),360*Nj);
% Hobs=reshape(Hobs,[3,Nj,360]);
% Recip1=sum(sum(Hobs.*ks,1),2);
% toc

[rt,kt]=resamplecoil(rs2,ks2,N,360);
[tri]=surftri(p',te2p');
[pp,~,tri2]=unique(tri(:));
tri=reshape(tri2,size(tri));
raux=zeros([size(rt),numel(pp)]);
kaux=zeros([size(kt),numel(pp)]);

for i=1:numel(pp)
    pt=p(:,pp(i));
    pt(3)=pt(3)+.005;
    nhat=pt/norm(pt);
[raux(:,:,i),kaux(:,:,:,i)]=movecoil(rt,kt,0,nhat,pt,360);
end
%%

tic
[Hobs]=computeHprimary(rv,jv,numel(rv)/3,reshape(raux,[3,numel(raux)/3]),numel(raux)/3);
Recip2=zeros([360 numel(pp)]);
%%
Hobs=reshape(Hobs,size(raux));
for j=1:numel(pp)
for i=1:360
Recip2(i,j)=sum(sum(Hobs(:,:,j).*kaux(:,:,i,j)));
end
end
toc
trisurf(tri,p(1,pp)',p(2,pp)',p(3,pp)',max(Recip2),'facecolor','interp','edgealpha',0);
xlabel('x')
zlabel('z')
ylabel('y')
axis equal
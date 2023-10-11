function modescript()
load ../example_point/siminibsetup0.mat
showfigs=1;
tri=surftri(p',te2p');
[ptri,~,tri]=unique(tri(:));
TR=triangulation(reshape(tri,[numel(tri)/3 3]),p(1,ptri)',p(2,ptri)',p(3,ptri)');
vn=vertexNormal(TR);
TR=triangulation(reshape(tri,[numel(tri)/3 3]),p(1,ptri)'+.004*vn(:,1),p(2,ptri)'+.004*vn(:,2),p(3,ptri)'+.004*vn(:,3));
clear ptri tri;
if showfigs==1
trisurf(TR,'edgealpha',0,'facealpha',.3,'facecolor',[141, 85, 36]/255);
light
lighting gouraud
hold on
%quiver3(TR.Points(:,1),TR.Points(:,2),TR.Points(:,3),vn(:,1),vn(:,2),vn(:,3));
end
robs=TR.Points';
[te2te]=gente2te(te2p);
[X,Y,Z]=ndgrid(-.005:.00025:.005,-.005:.00025:.005,0:-.0005:-.02);
r_roi=[X(:),Y(:),Z(:)];
[teid,bari]=pointlocation_c(r_roi',te2p-1,p,te2te,4);
teid=teid(teid~=0);
teid=unique(teid);
teid=teid(conductivity(teid)==0.2760);
if showfigs==1
tri=surftri(p',te2p(:,teid)');
trisurf(tri,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',1,'facecolor',[128, 128, 128]/255);
end

tecen=(p(:,te2p(1,teid))+p(:,te2p(2,teid))+p(:,te2p(3,teid))+p(:,te2p(4,teid)))/4;

%%%%%%%%begin ACA
Maxmodes=20;
capJ=ones([Maxmodes 1]);
norhist=zeros([Maxmodes 1]);
capI=ones([Maxmodes+1 1]);
v=zeros([numel(robs) Maxmodes]);
u=zeros([numel(tecen) Maxmodes]);
erri=zeros([Maxmodes 1]);
for i=1:Maxmodes
%compute row
[v(:,i),rv{i},jv{i}]=computerow(te2p,p,conductivity,capI(i),robs,teid);
if i~=1
    v(:,i)=v(:,i)-(u(capI(i),1:i-1)*v(:,1:i-1)')';
end
%define J
 C = setdiff(1:numel(robs),capJ);
pp=find(abs(v(C,i))==max(abs(v(C,i))));
capJ(i)=C(pp(1));
%compute column
    [u(:,i)]=computecolumn(te2p,p,conductivity,capJ(i),robs,tecen);
if i~=1
u(:,i)=u(:,i)-u(:,1:i-1)*v(capJ(i),1:i-1)';
end
    %how close are diagonal terms    
u(capI(i),i)
v(capJ(i),i)
v(:,i)=v(:,i)/v(capJ(i),i);
ukvk=norm(u(:,i))*norm(v(:,i));
norhist(i)=ukvk^2;
if i~=1
norhist(i)=norhist(i)+norhist(i-1)+2*sum(abs(u(:,1:i-1)'*u(:,i)).*abs(v(:,1:i-1)'*v(:,i)));
end
    %compute next index
    
 C = setdiff(1:numel(tecen),capI);
pp=find(abs(u(C,i))==max(abs(u(C,i))));
capI(i+1)=C(pp(1));
erri(i)=ukvk/sqrt(norhist(i));
er=ukvk/sqrt(norhist(i))
i
end

save modetest.mat;

capJ
end



function [Hobs,rv,jv]=computerow(te2p,p,conductivity,colid,robs,teid)
[dir,simid]=ind2sub([3 numel(teid)],colid)
that=zeros([3 1]);
that(dir)=1;
[rv,jv]=runcoderecip(te2p,p,conductivity,teid(simid)-1,that,2);
[Hobs]=computeHprimary(rv,jv,numel(rv)/3,robs,numel(robs)/3);
Hobs=-Hobs(:);
end


function [Efield]=computecolumn(te2p,p,conductivity,rowid,robs,tecen)
[dir,simid]=ind2sub([3 numel(robs)/3],rowid)
ks=zeros([3 1]);
ks(dir)=1;
[Efield]=runcodeks(te2p,p,conductivity,robs(:,simid),ks,tecen,2);
Efield=Efield(:);
end





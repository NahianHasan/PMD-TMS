
v1=zeros([3 Maxmodes+1]);

[dir,simid1]=ind2sub([3 numel(teid)],capI);
for i=1:Maxmodes+1
v1(dir(i),i)=1;
end

v2=zeros([3 Maxmodes]);
[dir,simid]=ind2sub([3 numel(robs)/3],capJ);
for i=1:Maxmodes
v2(dir(i),i)=1;
end

trisurf(TR,'edgealpha',0,'facealpha',.3,'facecolor',[141, 85, 36]/255);
light
lighting gouraud
hold on
axis equal
quiver3(robs(1,simid),robs(2,simid),robs(3,simid),v2(1,:),v2(2,:),v2(3,:));
hold on
quiver3(tecen(1,simid1),tecen(2,simid1),tecen(3,simid1),v1(1,:),v1(2,:),v1(3,:));


tri=surftri(p',te2p(:,teid)');
trisurf(tri,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',.3,'facecolor',[128, 128, 128]/255);
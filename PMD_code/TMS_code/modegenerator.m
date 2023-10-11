function modegenerator(simid,dire)

load /work/ljg24/loopstar/Arraycode/Simtools/reciprocitycodes/coilpositions.mat
%load ~/Downloads/coilpositions.mat

 v1=scalp.vertices(scalp.faces(:,2),:)-scalp.vertices(scalp.faces(:,1),:);
 v2=scalp.vertices(scalp.faces(:,3),:)-scalp.vertices(scalp.faces(:,1),:);
 nhat=cross(v1,v2,2);
 for i=1:numel(nhat(:,1))
 nhat(i,:)=nhat(i,:)/norm(nhat(i,:));
 end 
robs=(scalp.vertices(scalp.faces(:,1),:)...
     +scalp.vertices(scalp.faces(:,2),:)...
     +scalp.vertices(scalp.faces(:,3),:))/3+nhat*.005;
 
 
%determine volume currents
that=zeros([3,1]);
that(dire)=1;
[rv,jv]=runcoderecip(te2p,p,conductivity,teid(simid)-1,that,1);

%%


[Hobs]=computeHprimary(rv,jv,numel(rv)/3,robs',numel(robs)/3);
[Eobs]=computeEprimary(rv,jv,numel(rv)/3,robs',numel(robs)/3);

save(strcat('moderesults',num2str(simid),'_',num2str(dire),'.mat'),'Hobs','Eobs','rv','jv');

end
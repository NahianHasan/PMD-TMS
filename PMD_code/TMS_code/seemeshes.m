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
 A=zeros([numel(robs) 3*numel(teid)]);
 for dire=1:3
     for simid=1:numel(teid);
         
load(strcat('moderesults',num2str(simid),'_',num2str(dire),'.mat'),'Eobs');
A(:,dire+(simid-1)*3)=Eobs(:);
     end
 end
[u s v]=svd(A,0);
save solutionE.mat u s v A;

%%

te2p2=te2p(:,teid(:))';
for i=10
    col=reshape(v(:,i),[3,numel(teid)]);

col2=sqrt(sum(col.^2,1));
subplot(2,2,1),
tetramesh(te2p2,p',col2(:));
colorbar
title(strcat('Mode_{',num2str(i),'}'))
axis equal
view(2)
col3=reshape(uE(:,i),[3 numel(scalp.faces)/3]);
col4=sqrt(sum(col3.^2,1));
subplot(2,2,2),
trisurf(scalp.faces,scalp.vertices(:,1),scalp.vertices(:,2),scalp.vertices(:,3),col4(:),'edgealpha',0)
title('E-field')
axis equal
axis off
view(2)
col3=reshape(uH(:,i),[3 numel(scalp.faces)/3]);
col4=sqrt(sum(col3.^2,1));
subplot(2,2,3),
trisurf(scalp.faces,scalp.vertices(:,1),scalp.vertices(:,2),scalp.vertices(:,3),col4(:),'edgealpha',0)
title('H-field')
axis equal
axis off
view(2)
saveas(gcf,strcat('mode',num2str(i),'.tif'),'tiffn')
close all
end

semilogy(diag(s)/s(1))
ylabel('Mode efficiency')
xlabel('Mode id')
saveas(gcf,'svdplot.tif','tiffn')
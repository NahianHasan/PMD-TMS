function [Cout,rv,jv]=computerow(te2p,p,conductivity,colid,robs,teid,that1,that2,nhat)
[dir,simid]=ind2sub([3 numel(teid)],colid)
that=zeros([3 1]);
that(dir)=1;
[rv,jv]=runcoderecip(te2p,p,conductivity,teid(simid)-1,that,1);
[Aobs,curlAobs]=computeAprimary(rv,jv,numel(rv)/3,robs,numel(robs)/3);
Cout=zeros([numel(robs)/3 5]);
Cout(:,1)=Aobs(1,:).*that1(1,:)+Aobs(2,:).*that1(2,:)+Aobs(3,:).*that1(3,:);
Cout(:,2)=Aobs(1,:).*that2(1,:)+Aobs(2,:).*that2(2,:)+Aobs(3,:).*that2(3,:);
Cout(:,3)=Aobs(1,:).*nhat(1,:)+Aobs(2,:).*nhat(2,:)+Aobs(3,:).*nhat(3,:);
Cout(:,4)=curlAobs(1,:).*that1(1,:)+curlAobs(2,:).*that1(2,:)+curlAobs(3,:).*that1(3,:);
Cout(:,5)=curlAobs(1,:).*that2(1,:)+curlAobs(2,:).*that2(2,:)+curlAobs(3,:).*that2(3,:);
Cout=Cout(:);
end




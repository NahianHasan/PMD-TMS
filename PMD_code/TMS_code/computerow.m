function [Hobs,rv,jv]=computerow(te2p,p,conductivity,colid,robs,teid)
[dir,simid]=ind2sub([3 numel(teid)],colid)
that=zeros([3 1]);
that(dir)=1;
[rv,jv]=runcoderecip(te2p,p,conductivity,teid(simid)-1,that,1);
[Hobs]=computeHprimary(rv,jv,numel(rv)/3,robs,numel(robs)/3);
Hobs=-Hobs(:);
end




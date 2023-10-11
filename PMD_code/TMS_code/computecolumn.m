function [Efield]=computecolumn(te2p,p,conductivity,rowid,robs,tecen)
[dir,simid]=ind2sub([3 numel(robs)/3],rowid)
ks=zeros([3 1]);
ks(dir)=1;
[Efield]=runcodeks(te2p,p,conductivity,robs(:,simid),ks,tecen,1);
Efield=Efield(:);
end




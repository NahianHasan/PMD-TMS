function [Efield]=computecolumn(te2p,p,conductivity,rowid,robs,tecen,ncthat1,ncthat2,normal)
[simid,dir]=ind2sub([numel(robs)/5 5],rowid)
if dir==1  
js=ncthat1(:,simid);
[Efield]=runcode(te2p,p,conductivity,robs(:,simid),js,tecen,1);
Efield=Efield(:);
elseif dir==2
js=ncthat2(:,simid);
[Efield]=runcode(te2p,p,conductivity,robs(:,simid),js,tecen,1);
Efield=Efield(:);

elseif dir==3
elseif dir==4
ks=ncthat1(:,simid);
[Efield]=runcodeks(te2p,p,conductivity,robs(:,simid),ks,tecen,1);
Efield=Efield(:);

elseif dir==5
ks=ncthat2(:,simid);
[Efield]=runcode(te2p,p,conductivity,robs(:,simid),ks,tecen,1);
Efield=Efield(:);
end



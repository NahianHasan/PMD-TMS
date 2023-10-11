function [t2p,p]=figure8coilmodel(Nloop,ID,OD,wirehe,wireth,mesh_len)

% addpath('./wirecode')
NN=1000;
rs(:,1)=linspace(ID/2+wireth/2, OD/2-wireth/2, NN*Nloop).*cos(2*pi*(0:(Nloop*NN-1))/NN)-OD/2+wireth/2;
rs(:,2)=linspace(ID/2+wireth/2, OD/2-wireth/2, NN*Nloop).*sin(2*pi*(0:(Nloop*NN-1))/NN);
rs(:,3)=0;
rs(end+1:2*end,1)=-linspace(OD/2-wireth/2, ID/2+wireth/2, NN*Nloop).*cos(2*pi*(0:(Nloop*NN-1))/NN)+OD/2-wireth/2;
rs(NN*Nloop+1:end,2)=linspace(OD/2-wireth/2, ID/2+wireth/2, NN*Nloop).*sin(2*pi*(0:(Nloop*NN-1))/NN);
rs(NN*Nloop+1:end,3)=0;
rs(end-100:end,3)=linspace(0,wirehe+.00019,101);
 m=numel(rs)/3;
 rs(end+1:end+1000,1)=linspace(rs(end,1),rs(1,1),1000);
  rs(m+1:m+1000,2)=rs(m,2);
  rs(m+1:m+1000,3)=wirehe+.0002;

len=sqrt((rs(2:end,1)-rs(1:end-1,1)).^2+...
         (rs(2:end,2)-rs(1:end-1,2)).^2+...
         (rs(2:end,3)-rs(1:end-1,3)).^2);
         len=cumsum([0;len(:)]);
    N=ceil(len(end)/mesh_len);
    locs=0:len(end)/N:len(end);
    
     xcoil2=interp1(len,rs(1:end,1),locs,'line');
     ycoil2=interp1(len,rs(1:end,2),locs,'line');
     zcoil2=interp1(len,rs(1:end,3),locs,'line');

     %
% axis equal
clear rs
 rs(:,1)=xcoil2;
 rs(:,2)=ycoil2;
 rs(:,3)=zcoil2;

rs1=(rs(2:end,:)+rs(1:end-1,:))/2;
rs2=(rs(2:end,:)+rs(1:end-1,:))/2;

drs=rs(2:end,:)-rs(1:end-1,:);
%tan=drs(2:end,:)-drs(1:end-1,:);
%tan(end+1,:)=tan(end,:);
nhat=zeros(size(drs));
nhat(:,3)=1;
tan=cross(drs,nhat);
nh=1./sqrt(sum(drs.^2,2));
drs(:,1)=drs(:,1).*nh;
drs(:,2)=drs(:,2).*nh;
drs(:,3)=drs(:,3).*nh;
tan(:,3)=0;
nh=sqrt(sum(tan.^2,2));
nh=sign(nhat(:,3))./sqrt(sum(tan.^2,2));
tan(:,1)=tan(:,1).*nh;
tan(:,2)=tan(:,2).*nh;
tan(:,3)=tan(:,3).*nh;
nh=sign(nhat(:,3))./sqrt(sum(nhat.^2,2));
nhat(:,1)=nhat(:,1).*nh;
nhat(:,2)=nhat(:,2).*nh;
nhat(:,3)=nhat(:,3).*nh;

rs=rs1+wirehe*nhat/2;
ct=0;


rs1(:,1)=rs1(:,1)-wireth.*tan(:,1)/2;
rs1(:,2)=rs1(:,2)-wireth.*tan(:,2)/2;
rs1(:,3)=rs1(:,3)-wireth.*tan(:,3)/2;
rs2(:,1)=rs2(:,1)+wireth.*tan(:,1)/2;
rs2(:,2)=rs2(:,2)+wireth.*tan(:,2)/2;
rs2(:,3)=rs2(:,3)+wireth.*tan(:,3)/2;
rs3=rs1;rs4=rs2;
rs3=rs3+wirehe*nhat;
rs4=rs4+wirehe*nhat;

%plot3(rs1(:,1),rs1(:,2),rs1(:,3),'b')
%hold on
%plot3(rs2(:,1),rs2(:,2),rs2(:,3),'r')
%plot3(rs3(:,1),rs3(:,2),rs3(:,3),'g')
%plot3(rs4(:,1),rs4(:,2),rs4(:,3),'black')

[te2p,p,reg,port,portid,hexmesh]=maketetramesh(rs1,rs2,rs3,rs4);
% hold on
t2p=surftri(p,te2p);
% trisurf(t2p,p(:,1),p(:,2),p(:,3),'facecolor',[127 127 127]/256,'edgealpha',0,'facealpha',1);


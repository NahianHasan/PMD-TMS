function [te2p,p,reg,port,portid,hexmesh]=maketetramesh(rs1,rs2,rs3,rs4)
nl=numel(rs1(:,1));
p=cat(1,rs1,rs2,rs3,rs4);
hexmesh=zeros([nl-1 8]);
hexmesh(:,1)=1:nl-1;
hexmesh(:,2)=nl+1:2*nl-1;
hexmesh(:,3)=2*nl+1:3*nl-1;
hexmesh(:,4)=3*nl+1:4*nl-1;
hexmesh(:,5)=2:nl;
hexmesh(:,6)=nl+2:2*nl;
hexmesh(:,7)=2*nl+2:3*nl;
hexmesh(:,8)=3*nl+2:4*nl;

te2p=zeros([6*(nl-1) 4]);
%define te1
te2p(1:nl-1,1:4)=hexmesh(:,[5     6     7     3]);
te2p(nl:2*nl-2,1:4)=hexmesh(:,[7     6     8     3]);
te2p(2*nl-1:3*nl-3,1:4)=hexmesh(:,[2     6     1     3]);
te2p(3*nl-2:4*nl-4,1:4)=hexmesh(:,[5     1     6     3]);
te2p(4*nl-3:5*nl-5,1:4)=hexmesh(:,[3     6     8     4]);
te2p(5*nl-4:6*nl-6,1:4)=hexmesh(:,[2     6     3     4]);
reg=ones([6*(nl-1) 1]);
port=zeros([4 3]);
portid=[1;1;2;2];
port(1,:)=hexmesh(1,[1,2,3]);
port(2,:)=hexmesh(1,[2,4,3]);
port(3,:)=hexmesh(end,[6,7,8]);
port(4,:)=hexmesh(end,[5,6,7]);



end
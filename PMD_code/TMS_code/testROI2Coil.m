clc; clear all;
close all;
if ~isempty(dir('ACA_ADM.log'))
    diary off;
    delete('ACA_ADM.log');
end
diary ACA_ADM.log

addpath('Coil2ROI');

%%example setup parameters
mu0=1.25663706*10^-6;
eps0=8.85418782*10^-12;
dadt=6.664324407237550e+07;%coil current 
N=[17,17,2]; % the number of dipoles along x,y and z 17 by 17 by 2 is recommended (see publication
%for futher guidance)
FEMord=1;% order of the numerical approximation suggested 1 or 2 
th_hair=.005;%hair thickness or distance of coil from scalp
scth=.02;%width of the square desired scalp coil position search space
for ix=1:1
    name=['thirdordermesh',num2str(ix),'.mat'];
load(name); % thirdordermesh3.mat p te2p reg;
% load thirdordermesh0.mat p te2p conductivity;
te2p=te2p(1:4,:);
conductivity=.276*ones(size(reg));
[p2,~,te2p2]=unique(te2p(:));
p=p(:,p2);
te2p=reshape(te2p2,size(te2p));
p=p/max(p(:))*.085;

%load input data structures and make turn conducitivity into a tensor conductivity
%generate coil
subplot(1,2,1),
[rs ks]=genfig8(.056/2,.087/2,.006,9);%magstim specs
axis equal
rs=rs';ks=ks';
%turn conductivity into conductivity tensor
[~,conductivity]=ndgrid(1:9,conductivity(:));
conductivity([2,3,4,6,7,8],:)=0;
%generate ROI
roictr=[0 0 .07];
roirad=.005;
xx(conductivity(1,:)==.1260)=1; % White
% xx(conductivity(1,:)==.2760)=1; % Grey
[x y z]=ndgrid(0:pi/50:2*pi,0:pi/50:pi,0:roirad/50:roirad);
xp=z(:).*cos(x(:)).*sin(y(:))+roictr(1);
yp=z(:).*sin(x(:)).*sin(y(:))+roictr(2);
zp=z(:).*cos(y(:))+roictr(3);
clear x y z
TR=triangulation(te2p',p');
teid=pointLocation(TR,xp,yp,zp);
teid=unique(teid);
teid=teid(conductivity(1,teid)==.276);
clear xp yp zp;
that=zeros([3 numel(teid)]);
that(2,:)=1;

tri=surftri(p',te2p');
ROItri=surftri(p',te2p(:,teid)');
tri_matters=surftri(p',te2p(:,xx(:)==1)');
tri_matters2=surftri(p',te2p(:,conductivity(1,:)==0.126)');


[pp,Anor,tri]=findSurf(te2p,p,teid,scth,th_hair);


rs2(4,:)=1;
ncoil=size(pp,1);
nroi=size(teid,1);
tic;

ROI_id=1;
ZI1=ROI2coil_1xyz(tri,pp,Anor,dir0,te2p,p,conductivity,teid(ROI_id),rs,ks,omega,scth,th_hair,N,FEMord);
end

function ADM_WholeHead_cluster(k)

addpath('../maclibs');
addpath('../wirecode');
addpath('MatCode');

mu0=1.25663706*10^-6;
eps0=8.85418782*10^-12;
dadt=6.664324407237550e+07;%coil current 
N=[17,17,2]; % the number of dipoles along x,y and z 17 by 17 by 2 is recommended (see publication
%for futher guidance)
FEMord=1;% order of the numerical approximation suggested 1 or 2 
th_hair=.005;%hair thickness or distance of coil from scalp
scth=.05;%width of the square desired scalp coil position search space
load exampletest.mat p te2p conductivity;
%load input data structures and make turn conducitivity into a tensor conductivity
%generate coil
% subplot(1,2,1),
[rs ks]=genfig8(.056/2,.087/2,.006,9);%magstim specs
[t2pcoil,pcoil]=figure8coilmodel(9,.056,.087,.006,.0018,.001);
% axis equal
rs=rs';ks=ks';
%turn conductivity into conductivity tensor
[~,conductivity]=ndgrid(1:9,conductivity(:));
conductivity([2,3,4,6,7,8],:)=0;
%generate ROI
roilist=[1 2.5 5 10 15 20 25 30 35 40 45 50]*1e-3; % m
roictr=[0 0 -.015];
xyz={'x','y','z'};
dir00=1;  iroi=1;
roirad=roilist(iroi);
xx(conductivity(1,:)==.1260)=1; % White
xx(conductivity(1,:)==.2760)=1; % Grey
[x y z]=ndgrid(0:pi/150:2*pi,0:pi/150:pi,0:roirad/150:roirad);
xp=z(:).*cos(x(:)).*sin(y(:))+roictr(1);
yp=z(:).*sin(x(:)).*sin(y(:))+roictr(2);
zp=z(:).*cos(y(:))+roictr(3);
clear x y z
TR=triangulation(te2p',p');
teid=pointLocation(TR,xp,yp,zp);
teid=teid(isnan(teid)==0);
teid=unique(teid);
teid=teid(conductivity(1,teid)==.276);
clear xp yp zp;
clear teid
teid=1:numel(te2p)/4;
teid=teid(conductivity(1,teid)==.276);

that=zeros([3 numel(teid)]);
that(2,:)=1;

tri=surftri(p',te2p');
[ROItri,node5]=surftri(p',te2p(:,teid)');
tri_matters=surftri(p',te2p(:,xx(:)==1)');
tri_matters2=surftri(p',te2p(:,conductivity(1,:)==0.126)');

[pp,Anor,tri]=findSurf(te2p,p,teid,scth,th_hair);
pcoil(:,4)=1;
% idx = 12; view([0,100]); lightangle(gca,40,40);
% pcoilhead=(Anor(1:3,:,idx)*(pcoil'))';
% trisurf(t2pcoil,pcoilhead(:,1),pcoilhead(:,2),pcoilhead(:,3),'facecolor',[184 115 51]/256,'edgealpha',0,'facealpha',1);

teidN=teid; AnorN=Anor; ppN=pp;

disp(['Computing ROI2Coil ... ']);
tic
ROI2coil_Mat=ROI2coil_Matrix_new(k,te2p,p,tri,ppN,AnorN,conductivity,teidN,that,rs,ks,dadt,scth,th_hair,N,FEMord);
toc
fname=['ROI_Mat_',num2str(k),'.mat'];
save(fname,'ROI2coil_Mat','-v7.3');

disp(['Completed.']);

function ROI2coil_Mat=ROI2coil_Matrix_new(k,te2p,p,tri,ppN,AnorN,conductivity,teidN,that,rs,ks,dadt,scth,th_hair,N,FEMord)
rs2(4,:)=1;
ROI2coil_Mat0=[];
nt=size(teid,1);

disp(['ROI: ', num2str(k),' of ',num2str(nt)]);

ZIk=ROI2coil(tri,pp,Anor,te2p,p,conductivity,teid(k),rs,ks,omega,scth,th_hair,N,FEMord);
ROI2coil_Mat0=[ROI2coil_Mat0;reshape(ZIk,[],3)];
ROI2coil_Mat=ROI2coil_Mat0.'/omega;



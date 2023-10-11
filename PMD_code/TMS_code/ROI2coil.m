function [E_r2c]=ROI2coil(tri,pp,Anor,te2p,p,conductivity,teid,rs,ks,dadt,scth,th_hair,N,FEMord);

[E_aux]=genrecipmapxyzks(tri,pp,Anor,te2p,p,conductivity,teid,rs,ks,dadt,scth,th_hair,N,FEMord);

E_r2c=reshape(E_aux,[],3);
% E_r2c=reshape(E_r2c0.',[],1);
% ncoil=size(pp,1);
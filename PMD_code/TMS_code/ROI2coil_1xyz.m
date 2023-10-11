function [E_r2c]=ROI2coil_1xyz(tri,pp,Anor,dir0,te2p,p,conductivity,teid,rs,ks,dadt,scth,th_hair,N,FEMord);

[E_aux]=genrecipmapxyzks_ROI2coil(tri,pp,Anor,dir0,te2p,p,conductivity,teid,rs,ks,dadt,scth,th_hair,N,FEMord);

E_r2c=reshape(E_aux/dadt,[],1);
% % E_r2c=reshape(E_r2c0.',[],1);
% ncoil=size(pp,1);
function ROI2coil_Mat=ROI2coil_Matrix(te2p,p,tri,pp,Anor,conductivity,teid,that,rs,ks,omega,scth,th_hair,N,FEMord)

rs2(4,:)=1;
ROI2coil_Mat0=[];
nt=size(teid,1);
for k=1:nt
    disp(['ROI: ', num2str(k),' of ',num2str(nt)]);
    ZIk=ROI2coil(tri,pp,Anor,te2p,p,conductivity,teid(k),rs,ks,omega,scth,th_hair,N,FEMord);
%     ROI2coil_Mat=[ROI2coil_Mat;ZIk;];
    ROI2coil_Mat0=[ROI2coil_Mat0;reshape(ZIk,3,[])];
end
ROI2coil_Mat=ROI2coil_Mat0.'/omega;
nt;
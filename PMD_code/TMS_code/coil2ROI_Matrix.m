function coil2ROI_Mat=coil2ROI_Matrix(te2p,p,tri,pp,Anor,conductivity,teid,that,rs,ks,omega,scth,th_hair,N,FEMord)

rs2(4,:)=1;
ncoil=size(pp,1);

nt=size(teid,1);
for coilposid=1:ncoil
    for ang=1:360
        disp(['Coil # ',num2str(coilposid),' Angle: ', num2str(ang)]);
        ZIk=coil2ROI(coilposid,ang,rs,ks,Anor,teid,te2p,p,conductivity,FEMord);
        coil2ROI_Mat0(ang,coilposid,:)=ZIk;
    end
end
coil2ROI_Mat=reshape(coil2ROI_Mat0,360*ncoil,[]);

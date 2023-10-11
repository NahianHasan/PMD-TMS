function [Rs,Ks]=generateauxcoil(w,Anor,rs,ks,N)
    np=numel(Anor)/16;
    Nj=numel(rs)/3;
    rs2=zeros([3 Nj 360]);
    ks2=zeros([3 Nj 360]);
    for i=1:360
        phi=i*pi/180;
        rs2(3,:,i)=rs(3,:);
        ks2(:,:,i)=ks;
        rs2(1,:,i)= rs(1,:)*cos(phi)+rs(2,:)*sin(phi);
        rs2(2,:,i)=-rs(1,:)*sin(phi)+rs(2,:)*cos(phi);
    end
    [raux,kaux]=resamplecoil(rs2,ks2,N,360);
    raux(4,:)=1;
    %%generate dipole location array
    Rs=zeros([size(raux) np]);
    for i=1:np
        Rs(:,:,i)=Anor(:,:,i)*raux;
    end
    Rs=reshape(Rs(1:3,:,:),[3 prod(N)*np]);

    %%%generate dipoles
    Ks=zeros(3,prod(N),np);

    for i=1:360
        for j=1:np
            Ks(1,:,j)=Ks(1,:,j)+w(i+(j-1)*360)*kaux(3,:,i).*Anor(1,3,j);
            Ks(2,:,j)=Ks(2,:,j)+w(i+(j-1)*360)*kaux(3,:,i).*Anor(2,3,j);
            Ks(3,:,j)=Ks(3,:,j)+w(i+(j-1)*360)*kaux(3,:,i).*Anor(3,3,j);
        end
    end

    Ks=reshape(Ks(1:3,:,:),[3 prod(N)*np]);
end
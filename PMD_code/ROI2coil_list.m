function Eaux=ROI2coil_list(qvec,pp,Anor,te2p,p,conductivity,teid,rs,ks,N,FEMord)
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
    clear rs2 ks2;
    np=numel(pp)/3;
    %
    raux(4,:)=1;
    robs=zeros([size(raux) np]);
    for i=1:np
        robs(:,:,i)=Anor(:,:,i)*raux;
    end
    robs=reshape(robs(1:3,:,:),[3 prod(N)*np]);
    E_aux=zeros([1,360*np]);

    pA=p(:,te2p(1,:)); pB=p(:,te2p(2,:)); 
    pC=p(:,te2p(3,:)); pD=p(:,te2p(4,:)); 
    avec=pB-pC; bvec=pD-pC; cvec=pA-pC;

    vol=abs(dot(cross(avec,bvec,1),cvec,1))/6;
    clear pA pB pC pD avec bvec cvec
    vol=reshape(vol(teid),[1 numel(teid)]);
    tvol=sum(vol);
    tv=tvol./vol; tv3=repmat(tv,3,1);
    that0=reshape(qvec',3,[]);
    that=that0.*tv3;
    % for dir=1:3
    %     that=zeros([3 numel(teid)]);
    %     that(dir,:)=qvec(1:136)*tvol./vol;
    % end

    [rv11,jv11]=ROI2Coil_field(that,te2p,p,conductivity,teid-1,FEMord);

    [Hobs]=computeHprimary(rv11,jv11,numel(rv11)/3,robs,numel(robs)/3);

    Hobs=reshape(Hobs,[3 prod(N),np]);

    for j=1:np
        for i=1:360
            Eaux(i+360*(j-1))=-sum(kaux(3,:,i).*(Anor(1:3,3,j)'*Hobs(:,:,j)));%note, the calculated field is not multiplied by the coil current dadt
        end
    end
end


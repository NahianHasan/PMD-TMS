function [ux,vx,Misfit0,I,J]=ACA_ADM_ErrUVR(te2p,p,tri,pp,Anor,conductivity,teid,that,rs,ks,omega,scth,th_hair,N,FEMord)

% S0=svds(ROI2coil_Mat,80);
% S00=sort(S0,'descend');
% S01=S00/S00(1);

% raux(4,:)=1;
rs2(4,:)=1;
ncoil=size(pp,1);
nroi=size(teid,2);
tic;
I(1)=1; Zt=[]; %R=zeros(nroi*3,ncoil*360);
dir0=mod(I(1),nroi);
if dir0==0
    dir0=nroi;
end
ZI1=ROI2coil_1xyz(tri,pp,Anor,dir0,te2p,p,conductivity,teid(I(1)),rs,ks,omega,scth,th_hair,N,FEMord);
ZI1(isnan(ZI1)) = 0;
R1=ZI1'; %Z(I(1),:);
maxJk=find(abs(R1)==max(abs(R1)));
J(1)=maxJk(1);
vx(1,:)=R1/R1(J(1));
coilposid=floor(mod(J(1),360)/360)+1;
angid=mod(mod(J(1),360),360);
ZJ1=coil2ROI(coilposid,angid,rs,ks,Anor,teid,te2p,p,conductivity,FEMord);
R2=ZJ1; %Z(:,J(1));
ux(:,1)=R2;
maxIk=find(abs(R2)==max(abs(R2)));
I(2)=maxIk(1);
if ismember(I(2),I(1))
    R02=R2;
    R02(I(2))=[];
    I0=find(abs(R2)==max(abs(R02)));
    I(2)=I0(1);
end
k=1;
ZtN1=norm(ux(:,k))*norm(vx(k,:));
% Zt(:,:)=ux(:,:)*vx(:,:);
Misfit(1)=1;
idx0=[];
norm_uv=norm(ux(:,k))*norm(vx(k,:));
Misfit0(k)=norm_uv/ZtN1;
disp(['k = ', num2str(k),' MisfitUV: ', num2str(Misfit0(k))]);

for k=2:nroi*3
    ZtN0=ZtN1;
    % tic;
    uv=zeros(size(vx(k-1,:)));
    for l=1:k-1
        uv=uv+ux(I(k),l)*vx(l,:);
    end
    Ikid=floor(I(k)/3)+1;
    dir0=mod(I(k),3);
    if dir0==0
        dir0=3;
    end

    ZIk=ROI2coil_1xyz(tri,pp,Anor,dir0,te2p,p,conductivity,teid(Ikid),rs,ks,omega,scth,th_hair,N,FEMord);
    R1=ZIk.'-uv;
    [absR1,idx]=sort(abs(R1),'descend');
    Misfitk = NaN; kanan=0;
    while isnan(Misfitk)
        kanan = kanan + 1;
        maxIk=find(abs(R1)==max(abs(R1)));
        Jk0=maxIk(1);
        MisErr=1e-3; kErr=0; knull=0;
        while MisErr<3e-1
            kErr=kErr+1;
            if kErr==1
                J(k)=Jk0;
                JkErr(kErr)=J(k);
            else            
                J(k)=idx(kErr);
                kJkErr=kErr;
                while ismember(J(k),JkErr)
                    kJkErr=kJkErr+1;
                    J(k)=idx(kJkErr);
                end
                kErr0=kJkErr;
                while ismember(J(k),J(1:k-1))
                    kErr0=kErr0+1;
                    J(k)=idx(kErr0);
                end
            end
            JkErr(kErr)=J(k);

            vx(k,:)=R1/R1(J(k));
            vu=zeros(size(ux(:,k-1)));
            for l=1:k-1
                vu=vu+vx(l,J(k))*ux(:,l);
            end

            coilposid=floor(J(k)/360)+1;
            angid=mod(J(k),360);
            ZJk=coil2ROI(coilposid,angid,rs,ks,Anor,teid,te2p,p,conductivity,FEMord);
            R2=ZJk-vu;
            ux(:,k)=R2;
            uvn=0;
            maxJk=find(abs(R2)==max(abs(R2)));
            I(k+1)=maxJk(1);
            if ismember(I(k+1),I(1:k))
                R02=R2;
                R02(I(k+1))=[];
                I0=find(abs(R02)==max(abs(R02)));
                if I(k+1)>I0(1)
                    I(k+1)=I0(1);
                else
                    I(k+1)=I0(1)+1;
                end
            end
            norm_uv=norm(ux(:,k))*norm(vx(k,:));
            for l=1:k-1
        %         uvn=uvn+sum(ux(:,l).*ux(:,k))*sum(vx(l,:).*vx(k,:));
                uvn=uvn+sum(ux(:,l).*ux(:,k))*sum(vx(:,l).*vx(:,k));
            end
            ZtN1=sqrt(ZtN0^2+2*uvn+norm_uv^2);
            
        %     Zt(:,:)=ux*vx;
            
            Misfit0(k)=norm_uv/ZtN1;
            MisErr=Misfit0(k)/Misfit0(k-1);
            disp(['kErr: ',num2str(kErr),' Jk = ',num2str(J(k)),' MisErr: ',num2str(MisErr)]);
            if isnan(MisErr)
                break
            end
            
        end
        Misfitk = Misfit0(k);
    end
        
    disp(['k = ', num2str(k),' MisfitUV: ', num2str(Misfit0(k))]);
    if (isnan(Misfit0(k)))
        Misfit0(k)=[];
        ux(:,k)=[];
        vx(k,:)=[];
        break;
    end
    % time1=toc;
    clear uv vu
    if Misfit0(k)<1e-2
        break;
    end    
end
% ACA_Matrix=Zt;
k;


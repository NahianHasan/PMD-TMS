function ACA_Matrix=ACA_ADM_Fnorm(te2p,p,tri,pp,Anor,conductivity,teid,that,rs,ks,omega,scth,th_hair,N,FEMord)

% S0=svds(ROI2coil_Mat,80);
% S00=sort(S0,'descend');
% S01=S00/S00(1);

% raux(4,:)=1;
rs2(4,:)=1;
ncoil=size(pp,1);
nroi=size(teid,1);
tic;
I(1)=1; Zt=[]; R=zeros(nroi*3,ncoil*360);
dir0=mod(I(1),nroi);
ZI1=ROI2coil_1xyz(tri,pp,Anor,dir0,te2p,p,conductivity,teid(I(1)),rs,ks,omega,scth,th_hair,N,FEMord);
R(I(1),:)=ZI1'; %Z(I(1),:);
J(1)=find(R(I(1),:)==max(R(I(1),:)));
vx(1,:)=R(I(1),:)/R(I(1),J(1));
coilposid=floor(mod(J(1),360)/360)+1;
angid=mod(mod(J(1),360),360);
ZJ1=coil2ROI(coilposid,angid,rs,ks,Anor,teid,te2p,p,conductivity,FEMord);
R(:,J(1))=ZJ1; %Z(:,J(1));
ux(:,1)=R(:,J(1));
I(2)=find(R(:,J(1))==max(R(:,J(1))));
if ismember(I(2),I(1))
    R0=R;
    R0(I(2),:)=[];
    I0=find(R0(:,J(1))==max(R0(:,J(1))));
    if I0<=I(2)
        I(2)=I0;
    else
        I(2)=I0+1;
    end
end
%Zt(:,:)=ux(:,:)*vx(:,:);
k=1;
Misfit(1)=1;
idx0=[];
norm_uv=norm(ux(:,k))*norm(vx(k,:));
Misfit0(k)=norm_uv/norm(ux(:))/norm(vx(:));
% Misfit1(k)=norm(ROI2coil_Mat-Zt')/norm(ROI2coil_Mat);
% Misfit20=abs(ROI2coil_Mat-Zt')./max(max(ROI2coil_Mat));
% Misfit2(k)=max(max(Misfit20));
disp(['k = ', num2str(k),' MisfitUV: ', num2str(Misfit0(k))]); %, ' MisfitReal: ', num2str(Misfit1(k)), ' MisfitMax: ', num2str(Misfit2(k))]);
% disp(['k = ', num2str(k),' MisfitUV: ', num2str(Misfit0(k)), ' MisfitReal: ', num2str(Misfit1(k)), ' MisfitMax: ', num2str(Misfit2(k))]);
figure(100);
idx0=[idx0,1];
semilogy(idx0(k),Misfit0(k),'k.','markersize',12); hold on; grid on;
% semilogy(idx0(k),Misfit1(k),'r.','markersize',12); hold on; grid on;
% semilogy(idx0(k),Misfit2(k),'g.','markersize',12); hold on; grid on;
% semilogy(idx0(k),S01(k),'b.','markersize',12); hold on; grid on;
title('Convergence Curve');
xlabel('Rank'); ylabel('L^2 Norm Error')
h=legend;
for k=2:nroi*3
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
    if k==12
        k;
    end
    ZIk=ROI2coil_1xyz(tri,pp,Anor,dir0,te2p,p,conductivity,teid(Ikid),rs,ks,omega,scth,th_hair,N,FEMord);
    R(I(k),:)=ZIk.'-uv;
    maxIk=find(abs(R(I(k),:))==max(abs(R(I(k),:))));
    J(k)=maxIk(1);
    if ismember(J(k),J(1:k-1))
        R0=R;
        R0(:,J(k))=[];
        J0=find(R0(I(k),:)==max(R0(I(k),:)));
        if J0(1)<=J(k)
            J(k)=J0(1);
        else
            J(k)=J0(1)+1;
        end
    end
    vx(k,:)=R(I(k),:)/R(I(k),J(k));
    vu=zeros(size(ux(:,k-1)));
    for l=1:k-1
        vu=vu+vx(l,J(k))*ux(:,l);
    end
    
    coilposid=floor(J(k)/360)+1;
    angid=mod(J(k),360);
    ZJk=coil2ROI(coilposid,angid,rs,ks,Anor,teid,te2p,p,conductivity,FEMord);
    R(:,J(k))=ZJk-vu;
    ux(:,k)=R(:,J(k));    
    %Zt(:,:)=ux*vx;
    maxJk=find(abs(R(:,J(k)))==max(abs(R(:,J(k)))));
    I(k+1)=maxJk(1);
    if ismember(I(k+1),I(1:k))
        R0=R;
        R0(I(k+1),:)=[];
        I0=find(R0(:,J(k))==max(R0(:,J(k))));
        if I0(1)<=I(k+1)
            I(k+1)=I0(1);
        else
            I(k+1)=I0(1)+1;
        end
    end
    norm_uv=norm(ux(:,k))*norm(vx(k,:));
    Misfit0(k)=norm_uv/norm(ux(:))/norm(vx(:));
    %Misfit0(k)=norm_uv/norm(Zt);
%     Misfit1(k)=norm(ROI2coil_Mat-Zt')/norm(ROI2coil_Mat);
%     Misfit20=abs(ROI2coil_Mat-Zt')./max(max(ROI2coil_Mat));
%     Misfit2(k)=max(max(Misfit20));
    % time1=toc;
%     disp(['k = ', num2str(k),' MisfitUV: ', num2str(Misfit0(k)), ' MisfitReal: ', num2str(Misfit1(k)), ' MisfitMax: ', num2str(Misfit2(k))]);
    disp(['k = ', num2str(k),' MisfitUV: ', num2str(Misfit0(k))]); %, ' MisfitReal: ', num2str(Misfit1(k)), ' MisfitMax: ', num2str(Misfit2(k))]);
    idx0=[idx0,k];
    semilogy(idx0(k-1:k),Misfit0(k-1:k),'k-','linewidth',1.2); hold on;
    semilogy(idx0(k),Misfit0(k),'k.','markersize',12); hold on;    
%     semilogy(idx0(k-1:k),Misfit1(k-1:k),'r-','linewidth',1.2); hold on;
%     semilogy(idx0(k),Misfit1(k),'r.','markersize',12); hold on;    
%     semilogy(idx0(k-1:k),Misfit2(k-1:k),'g-','linewidth',1.2); hold on;
%     semilogy(idx0(k),Misfit2(k),'g.','markersize',12); hold on;    
%     semilogy(idx0(k-1:k),S01(k-1:k),'b-','linewidth',1.2); hold on;
%     semilogy(idx0(k),S01(k),'b.','markersize',12); hold on; 
    h.String={'ACA Stop Criterion'}; %,'ACA Error','SVD Error'};
    set(gca,'fontname','Times New Roman','fontsize',12);
    saveas(gcf,'ConvergenceCurve.png');
    saveas(gcf,'ConvergenceCurve.fig');
    % disp(['k = ', num2str(k)]);
    clear uv vu
    if Misfit0(k)<2e-2
        Zt(:,:)=ux*vx;
        break;
    end    
end
ACA_Matrix=Zt;
k;


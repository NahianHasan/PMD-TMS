function [Q,B,sapprox,Emap,sgt,PMDRES_2,PMDRES_F,RMSRE] = pmd_spherical_head(pmd_code_path,NModes,FEMORD,output_folder,sphere_density,subject_folder)
    addpath(fullfile(pmd_code_path,'TMS_code'));

    if ~exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]),'dir')
        mkdir(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]))
    end
    save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']);
    
    %Load a standard sphere and define the layer_source of spherical head
    [vv,ff] = icosphere(sphere_density);
    pcen=(vv(ff(:,1),:)+vv(ff(:,2),:)+vv(ff(:,3),:))/3;
    th=acos(pcen(:,3));
    ph=atan2(pcen(:,2),pcen(:,1));
    ff=ff(th<=(60/180*pi),:);%what part of the sphere to keep
    r_source=0.9*vv(vv(:,3)>1/2,:);
    k_source=vv(vv(:,3)>1/2,:);
    r_obs=0.7*vv;
    Ncoil=numel(r_source(:,1));
    Ntarg=numel(r_obs);
    r_source_tet = delaunay(r_source);
    r_source_tri = surftri(r_source,r_source_tet);
    
    %define dipoles
    coils.QP=r_source;
    coils.QN=k_source;
    coils.QW=ones(size(r_source(:,1)));
    coils.QPinds=(1:Ncoil)';
    coils.QPinds(:,2)=(1:Ncoil)';
    
    %analytical Efield
    Emap=zeros([Ntarg,Ncoil]);
    LFMm = tmslib_build_LFM_sphere([0,0,0], coils, r_obs);
    Emap=LFMm{1}';
    [~,sgt,~] = svd(Emap);
    sgt = diag(sgt);
    sgt = sgt(1:NModes);
    %PMD over the sphere
    g = gpuDevice();
    if numel(g) > 0
        Emap = gpuArray(Emap);
    end
    [Q,~]=qr(Emap*randn([Ncoil,NModes]),0);
    B=Q'*Emap;
    E_pred = Q*B;
    [~,sapprox,~]=svd(B,"econ");
    sapprox = diag(sapprox);
    if numel(g) > 0
        Q = gather(Q);
        B = gather(B);
        sapprox = gather(sapprox);
        E_pred = gather(E_pred);
        Emap = gather(Emap);
    end
    
    PMDRES_2=sapprox(NModes-4)/sapprox(1)*(1+sqrt((NModes-5)/(4))+exp(1)*sqrt(NModes)/5);
    PMDRES_F=sqrt(sum(sapprox(NModes-4).^2))/sqrt(sum(sapprox(1:NModes-4).^2))*(sqrt(1+(NModes-5)/4));
    RMSRE=vecnorm(E_pred-Emap,1)./vecnorm(Emap,2,1);
    save(save_file,'r_source','k_source','r_obs','r_source_tri','-v7.3')
    for ix=1:NModes
        Qi = Q(:,ix);
        Bi = B(ix,:);
        save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat']),'Qi','-v7.3');
        save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['B_',num2str(ix),'.mat']),'Bi','-v7.3');
    end
    if ~exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],'GT_E_Fields'),'dir')
        mkdir(fullfile(output_folder,['FEM_',num2str(FEMORD)],'GT_E_Fields'));
    end
    save(fullfile(output_folder,['FEM_',num2str(FEMORD)],'GT_E_Fields','E_org.mat'),'Emap','sgt','-v7.3')
    save(fullfile(output_folder,['FEM_',num2str(FEMORD)],'Errors_spherical_head.mat'),'RMSRE','PMDRES_2','PMDRES_F','-v7.3');
end
function [] = offline_parallel_stage_2(NModes,FEMORD,mode,output_folder)
    addpath(fullfile('.','TMS_code'));
    pathparts = strsplit(output_folder,filesep);
    subject_folder = pathparts{end};
    
    save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(mode),'.mat']);
    if exist(save_file,'file')
        variableInfo = who('-file', save_file);
        status=ismember('Efield', variableInfo);
    else
        if exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(mode),'.mat']),'file')
            variableInfo = who('-file', fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(mode),'.mat']));
            status=ismember('Qi', variableInfo);
        else
            status = 0;
        end
    end
    if ~status
        start_time = tic();
        %Load subject data
        F = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']));
        disp(['Coil2ROI = ',num2str(mode),'/',num2str(NModes)]);
        w = F.w;
        rcoil = F.rcoil;
        kcoil = F.kcoil;
        Anor_standardized = F.Anor_standardized;
        teid = F.teid;
        te2p = F.te2p;
        p = F.p;
        conductivity = F.conductivity;
        N_coil_dipoles = F.N_coil_dipoles;
        FEMORD = F.FEMORD;

        ro=(p(:,te2p(1,:))+p(:,te2p(2,:))+p(:,te2p(3,:))+p(:,te2p(4,:)))/4;
        [rsa,ksa]=generateauxcoil(w(:,mode),Anor_standardized,rcoil,kcoil,N_coil_dipoles);
        clear rcoil kcoil Anor_standardized N_coil_dipoles
        [te2p2,p2]=femgenmesh_c(te2p,p,FEMORD);
        [te2te]=gente2te(te2p);
        clear te2p p;
        %Step 2 assemble FEM matrix
        A=femassembleTensor(te2p2,p2,conductivity,FEMORD);
        %Step 3 generate right hand side of equation
        [rhs]=femgenrhskstensor(te2p2,p2,conductivity,rsa,ksa,FEMORD);
        clear conductivity
        % Step 4 delete one unknown and equation and define preconditioner
        A=A(1:end-1,1:end-1);
        PRECON=sparse(1:numel(A(:,1)),1:numel(A(:,1)),sqrt(full(diag(A))));
        rhs=rhs(1:end-1);
        % Step 5 solve system of equations
        [x,~,~,~,~,~]=minres(A,rhs,10^-10,10000,PRECON,PRECON);
        x(end+1)=0;
        clear A rhs PRECON
        % Step 6 evaluate field at desired locations
        Efield=FEMinterpolatorks(te2te,te2p2,p2,rsa,ksa,x,ro,FEMORD);
        Efield=Efield(:,teid);
        Efield=reshape(Efield,1,[]);
        Efield = Efield.';
        stage_2_Time = toc(start_time);
        if ~exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]),'dir')
            mkdir(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]))
        end
        save(save_file,'Efield','stage_2_Time');
    else
        disp(['Mode (Qi) = ',num2str(mode),' is already calculated'])
        return
    end
end


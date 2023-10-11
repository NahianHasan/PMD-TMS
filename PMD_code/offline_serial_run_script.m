function [] = offline_serial_run_script(pmd_code_path,NModes,msh_file,msh_file_read_fcn,msh_file_read_fcn_location,m2m_dir,...
                    FEMORD,output_folder,coil_model_file,th_hair,mapping_region,patch_angle,sphere_density,eeg_mni_source_file)
    start_time = tic;
    addpath(msh_file_read_fcn_location);
    addpath(fullfile(pmd_code_path,'TMS_code'));
    pathparts = strsplit(output_folder,filesep);
    subject_folder = pathparts{end};

    if ~exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]),'dir')
        mkdir(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]))
    end
    save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']);
    if ~exist(save_file,'file')
        [p,te2p,conductivity,reg,M] = load_msh_data(msh_file,msh_file_read_fcn);
        warning('off');
        [pp_standardized,Anor_standardized,tri_standardized,nhat_standardized] = generate_sample_coil_placement(msh_file,m2m_dir,...
                                                                                th_hair,patch_angle,sphere_density,eeg_mni_source_file);
        warning('on');
        %load the coil model
        [rcoil,kcoil,coil_model,coil_tri] = load_coil_model(pmd_code_path,coil_model_file,1);
        N_coil_dipoles=[17,17,2]; % the number of dipoles along x,y and z = 17 by 17 by 2 is recommended

        %select only those tetrahedrons within specified mapping region
        teid=1:numel(te2p)/4;
        if strcmpi(mapping_region,'GM')
            teid=teid(conductivity(1,teid)==.275);
        elseif strcmpi(mapping_region,'WM')
            teid=teid(conductivity(1,teid)==.126);
        elseif strcmpi(mapping_region,'GM+WM')
            teid=teid(conductivity(1,teid)==.275 | conductivity(1,teid)==.126);
        elseif strcmpi(mapping_region,'HEAD')
            teid = teid;
        else
            try
                teid = M.teid;
            catch
                disp("Please properly define the mapping region:'GM', 'WM', 'GM+WM', 'Head' or define the field teid inside the msh file.")
                return
            end
        end

        nc=size(pp_standardized,1)*360;
        w=randn(nc,NModes);%random matrix
        stage_1_Time = toc(start_time);
        save(save_file,'stage_1_Time','p','te2p','reg','conductivity','teid','NModes','FEMORD','pp_standardized','Anor_standardized',...
            'tri_standardized','w','nhat_standardized','rcoil','kcoil','coil_tri','coil_model','N_coil_dipoles','-v7.3');
    else
        load(save_file);
    end
    
    %%%%%%%%%%%%%%%%%%  Stage 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ix=1:NModes
        save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat']);
        if exist(save_file,'file')
            variableInfo = who('-file', save_file);
            status=ismember('Qi', variableInfo);
        else
            status = 0;
        end
        if ~status
            if ix>1
                disp('It seems some of the modes were deleted by someone unexpectedly. To correctly calculate the modes, starting again.')
            end
            for jx=1:NModes
                del_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['B_',num2str(jx),'.mat']);
                if exist(del_file,'file')
                    delete(del_file)
                end
            end
            break;
        end
    end
    if ~status
        Time_2e = 0;
        for ix=1:NModes
            save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(ix),'.mat']);
            start_time = tic();
            disp(['Coil2ROI = ',num2str(ix),'/',num2str(NModes)]);
            if exist(save_file,'file')
                variableInfo = who('-file', save_file);
                status=ismember('Efield', variableInfo);
            else
                status = 0;
            end
            if ~status
                ro=(p(:,te2p(1,:))+p(:,te2p(2,:))+p(:,te2p(3,:))+p(:,te2p(4,:)))/4;
                [rsa,ksa]=generateauxcoil(w(:,ix),Anor_standardized,rcoil,kcoil,N_coil_dipoles);
                [te2p2,p2]=femgenmesh_c(te2p,p,FEMORD);
                [te2te]=gente2te(te2p);
                %Step 2 assemble FEM matrix
                A=femassembleTensor(te2p2,p2,conductivity,FEMORD);
                %Step 3 generate right hand side of equation
                [rhs]=femgenrhskstensor(te2p2,p2,conductivity,rsa,ksa,FEMORD);
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
                Time_2 = toc(start_time);
                save(save_file,'Efield','Time_2');
                Q(:,ix) = Efield;
            else
                Q(:,ix) = load(save_file).Efield(:);
                Time_2 = load(save_file).Time_2;
            end
            Time_2e = Time_2e + Time_2;
        end
        stage_2_Time = Time_2e/NModes;

   
        %%%%%%%%%%%%%%%%%%  Stage 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]),'dir')
            mkdir(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]))
        end
        start_time = tic();
        %Take the QR decomposition of the mode matrix save the modes in the order of
        % highest to lowest singular values 
        [Q,~]=qr(Q,0);
        for ix=1:NModes
            Qi = Q(:,ix);
            save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat']),'Qi','-v7.3');
        end
        stage_3_Time = toc(start_time);
        save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']),'stage_2_Time','stage_3_Time','-append')
        %Delete the unnecessary mode files
        for ix=1:NModes
            del_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(ix),'.mat']);
            delete(del_file)
        end
    else
        for ix=1:NModes
            Q(:,ix) = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat'])).Qi;
        end
    end
    
    %%%%%%%%%%%%%%%%%%  Stage 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ix=1:NModes
        save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['B_',num2str(ix),'.mat']);
        if exist(save_file,'file')
            variableInfo = who('-file', save_file);
            status=ismember('Bi', variableInfo);
        else
            status = 0;
        end
        if ~status
            disp(['ROI2COIL = ',num2str(ix),'/',num2str(NModes)]);
            start_time = tic();
            Bi=ROI2coil_list((Q(:,ix))',pp_standardized,Anor_standardized,te2p,p,conductivity,teid,rcoil,kcoil,N_coil_dipoles,FEMORD);
            stage_4_Time = toc(start_time);
            save(save_file,'Bi','stage_4_Time','-v7.3');
        else
            continue;
        end
    end
    disp([newline,'Mode generration complete.',newline])
end
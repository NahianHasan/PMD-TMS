function [] = offline_parallel_stage_1(pmd_code_path,NModes,msh_file,msh_file_read_fcn,msh_file_read_fcn_location,m2m_dir,...
                    FEMORD,output_folder,coil_model_file,th_hair,mapping_region,patch_angle,sphere_density,eeg_mni_source_file)
    start_time = tic;
    addpath(msh_file_read_fcn_location);
    addpath(fullfile('.','TMS_code'));
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
        warning('on')
        %load the coil model
        [rcoil,kcoil,coil_model,coil_tri] = load_coil_model(pmd_code_path,coil_model_file,1);
        N_coil_dipoles=[17,17,2]; % the number of dipoles along x,y and z = 17 by 17 by 2 is recommended
        
        %select only those tetrahedrons within specified mapping region
        teid=1:numel(te2p)/4;
        if strcmpi(mapping_region,'GM')
            teid=teid(conductivity(1,teid)==.276);
        elseif strcmpi(mapping_region,'WM')
            teid=teid(conductivity(1,teid)==.126);
        elseif strcmpi(mapping_region,'GM+WM')
            teid=teid(conductivity(1,teid)==.276 | conductivity(1,teid)==.126);
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
    end
end


function [] = offline_parallel_stage_4(NModes,FEMORD,mode,output_folder)
    addpath(fullfile('.','TMS_code'));
    pathparts = strsplit(output_folder,filesep);
    subject_folder = pathparts{end};
    
    save_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['B_',num2str(mode),'.mat']);
    if exist(save_file,'file')
        variableInfo = who('-file', save_file);
        status=ismember('Bi', variableInfo);
    else
        status = 0;
    end
    
    if ~status
        start_time = tic;
        %Load each Q mode
        F = matfile(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']));
        disp(['ROI2COIL = ',num2str(mode),'/',num2str(NModes)]);
        Q = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(mode),'.mat'])).Qi;
        te2p = F.te2p;
        p = F.p;
        pp_standardized = F.pp_standardized;
        Anor_standardized = F.Anor_standardized;
        conductivity = F.conductivity;
        teid = F.teid;
        rcoil = F.rcoil;
        kcoil = F.kcoil;
        N_coil_dipoles = F.N_coil_dipoles;
        FEMORD = F.FEMORD;
        Bi=ROI2coil_list(Q',pp_standardized,Anor_standardized,te2p,p,conductivity,teid,rcoil,kcoil,N_coil_dipoles,FEMORD);
        if ~exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]),'dir')
            mkdir(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)]))
        end
        stage_4_Time = toc(start_time);
        save(save_file,'Bi','stage_4_Time','-v7.3');
    else
        disp(['Mode (Bi) = ',num2str(mode),' is already calculated'])
        return
    end
end
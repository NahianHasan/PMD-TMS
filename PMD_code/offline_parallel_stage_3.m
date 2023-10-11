function [] = offline_parallel_stage_3(NModes,FEMORD,output_folder)
    addpath(fullfile('.','TMS_code'));
    pathparts = strsplit(output_folder,filesep);
    subject_folder = pathparts{end};
    start_time = tic;
    for ix=1:NModes
        if exist(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat']),'file')
            variableInfo = who('-file', save_file);
            status=ismember('Qi', variableInfo);
        else
            status = 0;
        end
        if status==0
            if ix>1
                disp('It seems some of the modes were deleted by someone unexpectedly.')
            end
            break;
        end
    end
    if ~status
        %Load the spanning modes
        Time_2 = 0;
        for i=1:NModes
            disp(['Reading Modess = ',num2str(i),'/',num2str(NModes)]);
            Q(:,i) = load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(i),'.mat'])).Efield(:);
            Time_2 = Time_2 + load(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(i),'.mat'])).stage_2_Time;
        end
        stage_2_Time = Time_2/NModes;
        %Take the QR decomposition of the mode matrix save the modes in the order of
        % highest to lowest singular values 
        [Q,~]=qr(Q,0);
        for ix=1:NModes
            Qi = Q(:,ix);
            save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['Q_',num2str(ix),'.mat']),'Qi','-v7.3');
        end
        %Delete the unnecessary mode files
        for i=1:NModes
            del_file = fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],['coil2roi_',num2str(i),'.mat']);
            delete(del_file)
        end
        stage_3_Time = toc(start_time);
        save(fullfile(output_folder,['FEM_',num2str(FEMORD)],['Modes_',num2str(NModes)],[subject_folder,'_FEM_',num2str(FEMORD),'.mat']),'stage_2_Time','stage_3_Time','-append')
    end
end

function [rcoil,kcoil,coil_model,coil_tri] = load_coil_model(pmd_code_path,coil_model_file,save_model)
    clear rcoil kcoil
    if strcmpi(coil_model_file,'fig-8')
        [rcoil,kcoil]=genfig8(.056/2,.087/2,.006,9);%magstim specs
        rcoil=rcoil';kcoil=kcoil';
        [coil_tri,~]=figure8coilmodel(9,.056,.087,.006,.0018,.001);
        coil_model = 'Figure-8-70-mm';
    else
        M = load(coil_model_file);
        t = strsplit(coil_model_file,filesep);
        coil_model = t{end}(1:end-4);
        rcoil = M.rcoil;
        kcoil = M.kcoil;
        try
            coil_tri = M.tri;
        catch
            coil_tri = [];
            disp('coil triangulation info is not provided')
        end
    end
    if save_model
        if ~exist(fullfile(pmd_code_path,'coil_data.mat'),'file')
            save(fullfile(pmd_code_path,'coil_data.mat'),'coil_model','rcoil','kcoil','coil_tri','-v7.3')
        end
    end
end
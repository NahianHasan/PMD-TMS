function [] = main_offline_stage(pmd_code_path,msh_file,msh_file_read_fcn,msh_file_read_fcn_location,m2m_dir,NModes,FEMORD,...
                                output_folder,run_mode,cluster_parameter_file,th_hair,options)
    arguments
        pmd_code_path char
        msh_file char
        msh_file_read_fcn char
        msh_file_read_fcn_location char
        m2m_dir char
        NModes double
        FEMORD double
        output_folder char
        run_mode char
        cluster_parameter_file char
        th_hair double
        options.mapping_region char = "GM+WM"
        options.patch_angle (1,1) {mustBeNumeric} = 60
        options.sphere_density (1,1) {mustBeNumeric} = 5
        options.eeg_mni_source_file char = "EEG10-10_UI_Jurak_2007.csv"
        options.coil_model_file char = "fig-8"
    end
    
    if strcmpi(run_mode,'parallel')
        [cluster_params] = collect_cluster_parameters(cluster_parameter_file);
        cluster_name = cluster_params{1};
        stage_1_cpu = cluster_params{2};
        stage_2_cpu = cluster_params{3};
        stage_3_cpu = cluster_params{4};
        stage_4_cpu = cluster_params{5};
        stage_1_max_walltime = cluster_params{6};
        stage_2_max_walltime = cluster_params{7};
        stage_3_max_walltime = cluster_params{8};
        stage_4_max_walltime = cluster_params{9};
        matlab_module_version = cluster_params{10};
    end
    if ~exist(output_folder,'dir')
        mkdir(output_folder)
    end
    if strcmpi(run_mode,'parallel')
        if ~exist(fullfile('.','slurm_output'),'dir')
            mkdir(fullfile('.','slurm_output'))
        end
        parallel_bash_script_name = fullfile(pmd_code_path,'offline_parallel_run_script.sh');
        if ispc
            conversion_command = 'unix2dos';
        else
            conversion_command = 'dos2unix';
        end
        cmdStr = [conversion_command ' ' parallel_bash_script_name];
        system(cmdStr);
        for stage=1:4
            cmdStr = [conversion_command,' ',fullfile(pmd_code_path,['offline_parallel_stage_',num2str(stage),'.sh'])];
            system(cmdStr);
        end
    
        cmdStr = ['bash' ' ' parallel_bash_script_name ' ' pmd_code_path ' ' num2str(NModes) ' ' msh_file ' '...
                msh_file_read_fcn ' ' msh_file_read_fcn_location ' ' m2m_dir ' ' num2str(FEMORD) ' ' output_folder ' '...
                options.coil_model_file ' ' num2str(th_hair) ' ' options.mapping_region ' ' num2str(options.patch_angle) ' ' ...
                num2str(options.sphere_density) ' ' options.eeg_mni_source_file ' ' stage_1_cpu ' ' stage_2_cpu ' '...
                stage_3_cpu ' ' stage_4_cpu ' ' cluster_name ' ' stage_1_max_walltime ' ' stage_2_max_walltime ' '...
                stage_3_max_walltime ' ' stage_4_max_walltime ' ' matlab_module_version];
        system(cmdStr);
        disp([newline,'Modes are being calculated in the cluster',newline])
    elseif strcmpi(run_mode,'serial')
        offline_serial_run_script(pmd_code_path,NModes,msh_file,msh_file_read_fcn,msh_file_read_fcn_location,m2m_dir,...
                    FEMORD,output_folder,options.coil_model_file,th_hair,options.mapping_region,...
                    options.patch_angle,options.sphere_density,options.eeg_mni_source_file);
    else
        disp("Please specify whether to run the offline stage in 'parallel' or 'serial'")
    end
end

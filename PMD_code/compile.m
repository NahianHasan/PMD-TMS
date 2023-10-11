function [output_folder] = compile(pmd_code_path,msh_file_read_fcn_location,output_directory,subject_folder,varargin)
    if ~isempty(varargin)
        simnibs_installation_directory = varargin{1};
        addpath(simnibs_installation_directory);
        addpath(fullfile(simnibs_installation_directory,'matlab'));
    end
    addpath(pmd_code_path)
    addpath(fullfile(pmd_code_path,'TMS_code'));
    addpath(msh_file_read_fcn_location);
    output_folder=fullfile(output_directory,subject_folder);
    
end
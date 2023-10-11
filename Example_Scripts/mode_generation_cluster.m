clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining Paths
pmd_code_path = '/scratch/bell/hasan34/data/Proj-2_PMD/Code_Github/PMD_code';
addpath(pmd_code_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter Set-up
NModes=110;%number of modes
msh_file='/scratch/bell/hasan34/data/head_model_data/mri2mesh_models/Updated/BA40/mri2mesh_models/S35/S35.msh';
msh_file_read_fcn = 'mesh_load_gmsh4';%useful for '.msh' msh-files. Also make sure the function is in the matlab search path. For '.mat' msh-files, it's not necessary.
msh_file_read_fcn_location = fullfile('/home/hasan34/SimNIBS-3.2/matlab');
m2m_dir = '/scratch/bell/hasan34/data/head_model_data/mri2mesh_models/Updated/BA40/mri2mesh_models/S35/m2m_S35';
FEMORD=1;%FEM order = 1,2, or, 3
output_directory = '/scratch/bell/hasan34/data/Proj-2_PMD/Code_Github/solutions_sample';%output directory for the results
subject_folder = 'BA40_S35_mri2mesh';%name of the subject
run_mode='parallel';%options = 'serial','parallel' (for HPC clusters);
%If parallel, provide the cluster parameters in a separate csv file (cluster_parameters.csv) (compatible with slurm scripting)
cluster_parameter_file = '/scratch/bell/hasan34/data/Proj-2_PMD/Code_Github/Example_Scripts/cluster_parameters.csv';
th_hair = 0.005;%distance of coil bottom center from scalp to accommodate hair thickness
simnibs_installation_directory = '/home/hasan34/SimNIBS-3.2';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mode Calculation
[output_folder] = compile(pmd_code_path,msh_file_read_fcn_location,output_directory,subject_folder,simnibs_installation_directory);
main_offline_stage(pmd_code_path,msh_file,msh_file_read_fcn,msh_file_read_fcn_location,m2m_dir,NModes,FEMORD,...
                        output_folder,run_mode,cluster_parameter_file,th_hair)
             
%% Optional name-value argument list to be provided to "main_offline_stage()"
%1. sphere_density = 5;%to specify the density of the standardized EEG points density. 
%                         The higher the number, the denser the points. default=5.
%2. patch_angle = 60;%The standardized EEG surface patch angle within which the points are sampled. default=60.
%
%3. eeg_mni_source_file = 'EEG10-10_UI_Jurak_2007.csv';%define the name of standard eeg-coordinate source file inside the
%                                                               m2m folder. Default is 'EEG10-10_UI_Jurak_2007.csv'. 
%                                                               If something different is required, provide the name of the file that
%                                                               is inside the ElectrodeCaps_MNI folder.
%4. coil_model_file = '';%specify the coil model file. If not provided, a Fig-8 coil model will be used.
%5. mapping_region = 'GM';%options = ['GM','WM','GM+WM','head',''];If nothing is specified, 
%                               provide a list of teids inside the msh file as a separate field; default = 'GM+WM';
                                                            
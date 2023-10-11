clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining Paths
pmd_code_path = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.1_PMD/Codes/Code_Github/PMD_code';
addpath(pmd_code_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter Set-up
NModes=110;%number of modes
FEMORD=1;%FEM order = 1,2, or, 3
sphere_density = 5;%to specify the density of the standardized EEG points density. 
%                         The higher the number, the denser the points. default=5.
msh_file_read_fcn = 'mesh_load_gmsh4';%useful for '.msh' msh-files. Also make sure the function is in the matlab search path. For '.mat' msh-files, it's not necessary.
msh_file_read_fcn_location = fullfile('/home/wang3007-1/SimNIBS-3.2/matlab');
output_directory = '/media/wang3007-1/Nahian_5TB/Research/Purdue_University/Professor_Gomez/Projects/Proj-2.1_PMD/Codes/Code_Github/solutions_sample';%output directory for the results
subject_folder = 'spherical_head';%name of the subject

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mode Calculation
[output_folder] = compile(pmd_code_path,msh_file_read_fcn_location,output_directory,subject_folder);
[Q,B,sapprox,Emap,sgt,PMDRES_2,PMDRES_F,RMSRE] = pmd_spherical_head(pmd_code_path,NModes,FEMORD,output_folder,sphere_density,subject_folder);
     
%%%%%%%%%%%%%%%
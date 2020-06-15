function [imag_vol, path, files] = DICOMLoad(start_path)
%% Reads in segmented lung masks that were exported as DICOM files from Amira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by I.R. Stecker 05/29/2020
% 
% Input:    start_path: initial path before changing directories (i.e., pwd)
%
% Output:   imag_vol:   3D image output volume
%           FileNames:  names of all files used to generate volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Checking nargin and setting up paths to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check number of input arguments
if nargin == 0
    start_path = pwd;
end

% Clear old data
if exist('files')~=0
    clear('files');
end

% Select folder containing exported DICOM files from Amira
path = uigetdir('Select Exported Dicom Folder from Amira');
files = dir(fullfile(path, '*.dcm'));
cd(path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparing to load in files to create segmented mask image from Amira
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the total number of slices found in folder
num_files = size(files,1);
num_slices= num_files;
msg = sprintf(['\n' int2str(num_slices) ' slices found.']);
disp(msg);

% Load image slices into a cell
all_slices = cell(num_files,1);

for indx = 1:num_files
    all_slices{indx} = dicomread(files(indx).name);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate 3D volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create mask volume output variable
imag_size = size(all_slices{1});
imag_vol = zeros(imag_size(1),imag_size(2),num_slices);

for slice = 1:num_slices
    imag_vol(:,:,slice) = all_slices{slice};
end

cd(start_path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
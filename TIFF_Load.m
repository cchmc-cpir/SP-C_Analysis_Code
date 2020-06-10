function [imag_vol, path, FileNames] = TIFF_Load(save_name);
%% A funtion of read in lung masks segmented in Amira and exported as TIF files
%
% Written by Z.I. Cleveland 03/26/2015
% 
% Input:    save_save:   command to save as a .mat file. 
%
% Output:   imag_vol:   3D image output volume
%           FileNames:  names of all files used to generate volume
%
%% 
%%

% Load data
if exist('FileNames')~=0 % clear memory if needed
    clear('FileNames');
end

start_path = pwd;
[FileNames,path]=uigetfile('*.tif','Select TIFF files','multiselect','on');
cd(path);

num_files=length(FileNames);
num_slices= num_files;
msg = sprintf(['\n' int2str(num_slices) ' slices found.']);
disp(msg);

% load image slices
all_slices = cell(num_files,1);

for indx = 1:num_files
    img_slice = imread(FileNames{indx});
    all_slices{indx} = img_slice;
end

% generate 3D volume
imag_size = size(all_slices{1});
imag_vol = zeros(imag_size(1),imag_size(2),num_slices);

for slice = 1:num_slices
    imag_vol(:,:,slice) = all_slices{slice};
end

% Save data
if nargin >=1
    
    save_check = ischar(save_name);
    
    if save_check == 1 
        save_name = [save_name '.mat'];
        cd ../        
        save(save_name);
        cd(path);
    else
        msg = sprintf('\nUnrecognized save command. Data was not saved.');
        disp(msg);
    end
end

cd(start_path)

end
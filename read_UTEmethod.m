function [MethodParams] = read_UTEmethod(path, zerofilling)
%% Radial UTE Reconstruction Protocol

% Original Author: Jinbang Guo; Co-Author: Ian Stecker
% This is a re-worked function of the ReadTextData.m file 

% The purpose of this function is to re-calculate the trajectories of the
% measured k-space data that was calculated by Bruker to be used later for 
% reconstruction for inspiration/expiration

%% Inputs of Function:
    % path: working folder containing data to retrospectively gate(typically, path = pwd)
    
    % zerofilling: must be string that contains either:
    % 'yes' or 'no'
        % If no, does not implement zero filling for reconstructed image (anisotropic)
        % If yes, implements a zero filling for reconstructed image(isotropic)
    
    % method_file: must be string that contains the name of the specific
    % method file to be loaded into the workspace (i.e. 'method', 'method_water')
    
    % acqp_file: must be string that contains the name of the specific acqp
    % file to be loaded into the workspace (i.e. 'acqp', 'acqp_retro')
    
    % traj_filename: must be string that contains the name of desired
    % filename you want to call the recalculated trajectory file
    % (i.e. 'traj_measured', 'traj_water')
    
    % savefile: must be string that contains the name of desired filename
    % you want to call the saved MATLAB data for filing purposes
    % (i.e. 'bleo03_readtraj', or whatever you wish to call it)

%% Outputs of Function:
    % The re-calculated trajectories of radial UTE where the file
    % is named after traj_filename. The file is saved in the current mouse folder
    
    % The MATLAB data that stores the path of the working directory of the
    % mouse folder as well as if the image was zero-filled or not which is
    % named after savefile. The file is saved in current mouse folder 

%% Sequence of Protocol: 
%   1. CalcTraj       - Creates the 'traj_measured' file
%
%   2. RetroGate      - Retrospecitively gates FID magnitude
%                       Generates 'fid_(insp/exp)'  file
%                       Generates 'traj_(insp/exp)' file
%
%   3. sdc_retro      - Loads in generated files; recons
%                       either 'img_(insp/exp).raw' file

%% Prints directory starting at:
start_path = pwd;

%% Selecting the Desired Method File from Given Folder
% Function to read Method File and get trajectories

% If path is passed through the function, use that path, otherwise have user
% select the method file from the trajectory measurement
if nargin == 0
    [filename,path]=uigetfile('*','Select Method file');
    %pfile_name = strcat(path, filename);
    cd(path)
    meas_method = filename;
    theo_method = 'method';
else
    cd(path)
    methodfiles = dir('method*');
    if numel(methodfiles.name)>6
        [method1,method2] = methodfiles.name;
        if length(method1) == 6
            theo_method = method1;
            meas_method = method2;
        else
            theo_method = method2;
            meas_method = method1;
        end
    else
        method = methodfiles.name;
        theo_method = method;
        meas_method = method;
    end
end
%% Obtaining Relevant Scan Parameters from Mouse Scan Method File
% Finding NPoints, ACQShift, NPro, NTE
scan_method = fopen(theo_method);
methodread = textscan(scan_method,'%s','delimiter','\n'); 
methodread = methodread{1};

for index = 1:size(methodread, 1)
        test_str = char(methodread{index});
        
        if length(test_str) > 17
            if strcmp(test_str(1:18), '##$PVM_TrajSamples') == 1
                readout_str = methodread{index};
                methodread_NPoints = readout_str(20:length(readout_str));
            end
        end
        
        if length(test_str) > 5
            if strcmp(test_str(1:6),'##$NTE') == 1
                readout_str = methodread{index};
                methodread_NTE = readout_str(8:length(readout_str));
            end
        end
        
        if length(test_str) > 6
            if strcmp(test_str(1:7), '##$NPro') == 1
                readout_str = methodread{index};
                methodread_NPro = readout_str(9:length(readout_str));
            end
        end
        
        if length(test_str) > 10
            if strcmp(test_str(1:11), '##$AcqShift') == 1
                readout_str = methodread{index}; 
                methodread_AcqShift = readout_str(13:length(readout_str));
            end
        end
        
        if length(test_str) > 15
            if strcmp(test_str(1:16), '##$PVM_EncMatrix') == 1
			readout_str = methodread{index+1}; 
			methodread_ImageMatrix = readout_str(1:length(readout_str));
            end
        end
end

fclose(scan_method);

% Returning desired parameters from method file
%NPointsarray = methodread_NPoints;
NPoints         = str2num(methodread_NPoints);
NShift          = str2num(methodread_AcqShift);
NTE             = str2num(methodread_NTE);
PreImageMatrix  = str2num(methodread_ImageMatrix);
NPro_tot        = str2num(methodread_NPro);

% Dividing by Total # of TEs
NPro = NPro_tot / NTE;

% Real Image Matrix Size
ImageMatrix = PreImageMatrix - 2*NShift;

% Interleaves = 13; % for keyhole trajectory
for k = 1:numel(scan_method)
    
    scan_method = fopen(meas_method);
    
    % 'goldmean' for Golden Mean Trajectory, 'keyhole' for Keyhole Trajectory
    trajmode = 'goldmean';
    
    % For 2D Golden Mean Trajectory
    phi = [0.46557123 0.68232780];
    
    % Obtaining the lines of forced trajectory measurements for Kx, Ky, Kz
    line = fgetl(scan_method);
    
    while(~strcmp(strtok(line,'='), '##$PVM_TrajKx'))
        line = fgetl(scan_method);
    end
    
    % Creating NPoints x 1 zeros matrix to be filled by traj from methods file
    TrajKx = zeros(NPoints, 1);
    KxCount = 0;
    while(~strcmp(strtok(line,'='), '##$PVM_TrajKy'))
        line = fgetl(scan_method);
        temp = (str2num(line))';
        TrajKx(KxCount + 1:size(temp, 1)+ KxCount) = temp;
        KxCount = KxCount + size(temp, 1);
    end
    
    TrajKy = zeros(NPoints,1);
    KyCount = 0;
    while(~strcmp(strtok(line, '='), '##$PVM_TrajKz'))
        line = fgetl(scan_method);
        temp = (str2num(line))';
        TrajKy(KyCount + 1:size(temp, 1) + KyCount) = temp;
        KyCount = KyCount + size(temp, 1);
    end
    
    TrajKz = zeros(NPoints, 1);
    KzCount = 0;
    while(~strcmp(strtok(line, '='), '##$PVM_TrajBx'))
        line = fgetl(scan_method);
        temp = (str2num(line))';
        TrajKz(KzCount + 1:size(temp, 1) + KzCount) = temp;
        KzCount = KzCount + size(temp, 1);
    end
    fclose(scan_method);
    
    maxkx = max(TrajKx);
    maxky = max(TrajKy);
    maxkz = max(TrajKz);
    
    if strcmp(trajmode, 'keyhole')
        nviews = NPro;
        keys = Interleaves;
        halfnviews = int32((nviews - 1) / 2);
        keyviews = nviews/keys;
        sf = int32((keyviews - 1) / 2);
        primeplus = 203;
        r = zeros(nviews, 1);
        p = zeros(nviews, 1);
        s = zeros(nviews, 1);
        fl = -1;
        grad_indx=0;
        for j = 0:(keys - 1)
            for i = 1:keyviews
                indx = j + (i - 1) * keys;
                f = 1 - double(indx) / double(halfnviews);
                if(f < -1)
                    f = -1;
                end
                ang = primeplus * indx * pi / 180;
                d = sqrt(1 - f * f);
                grad_indx = grad_indx + 1;
                r(grad_indx) = d * cos(ang);
                p(grad_indx) = d * sin(ang);
                if(i <= sf)
                    s(grad_indx) = sqrt(1 - d * d);
                else
                    s(grad_indx) = fl * sqrt(1 - d * d);
                end
            end
        end
    else
        nviews = NPro;

        r = zeros(nviews, 1);
        p = zeros(nviews, 1);
        s = zeros(nviews, 1);
        for i = 1:nviews
            s(i) = 2 * mod((i - 1) * phi(1), 1) - 1;
            alpha = 2 * pi * mod((i - 1) * phi(2), 1);
            d=sqrt(1-s(i)^2);
            r(i) = d * cos(alpha);
            p(i) = d * sin(alpha);
        end
    end
    %% Updated by Ian 10/17/17: Created setting to do trajectories for
    %% an/isotropic resolution (Zero Filling)
    if  strcmp(zerofilling, 'no')
        %These 3 lines indicate anisotropic, i.e. no zero filling
        TrajKx = TrajKx / maxkx / 2;
        TrajKy = TrajKy / maxky / 2;
        TrajKz = TrajKz / maxkz / 2;
        
        %These lines indicate anisotropic
        trajectory = zeros(3, NPoints, nviews);
        for i = 1:nviews
            trajectory(1, :, i) = r(i) * TrajKx;
            trajectory(2, : ,i) = p(i) * TrajKy;
            trajectory(3, :, i) = s(i) * TrajKz;
        end
    else
        %These lines indicate isotropic
        trajectory = zeros(3, NPoints, nviews);
        for i = 1:nviews
            trajectory(1, :, i) = r(i) * TrajKx;
            trajectory(2, :, i) = p(i) * TrajKy;
            trajectory(3, :, i) = s(i) * TrajKz;
        end
    end
    
%% IS Note 4/13/2020 
% % You can add an extra I/P variable called traj_filename, just give a
% % string name and it will write a trajectory file to the folder

%     % writes the re-calculated trajectories as its own file
%     scan_method = fopen(fullfile(path, traj_filename),'w');
%     fwrite(scan_method, trajectory, 'double');
%     fclose(scan_method);
end

%% Optional: Can Save the Workspace as a MATLAB Data File
MethodParams.ImageMatrix = ImageMatrix; 
MethodParams.NPoints = NPoints;
MethodParams.NPro    = NPro;
MethodParams.NShift  = NShift;
MethodParams.NTE     = NTE;
MethodParams.Traj = trajectory;
MethodParams.ZeroFilling = zerofilling;

% save(savefile, 'MethodParams')
cd(start_path)
end


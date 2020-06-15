function [ImgResults] = ImgAnalysis_SPC(path, roi_val, lung_val, savefile)
%% ImgAnalysis_SPC performs all computations that yielded results in manuscript 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by I.R. Stecker 05/29/2020

% Purpose of this function: to calculate the mean and cv of the whole
% lung volume. In addition, it will also measure the percentage of high 
% signal volume from segmented lung masks and also the tidal volume
% from the selected folder. The folder should contain both the 
% reconstructed image and file(s) from manual/automatic segmentations of
% the expiratory and inspriatory lung volumes

% Note: There are multiple ways to load in both the image & segmented
%       masks. Here, I primarily focus on having my images stored as TIFF
%       and my segmented masks stored as DICOM.

%       Alternative ways to store images can be by: .mat, .raw, .am, .nii
%       Alternative ways to store segmented masks can be by: .am, .nii

% Input:        path: initial path before changing directories (i.e., pwd)

%               roi_val: this is the assigned mask label from
%               Amira/ITK_Snap segmentations. This mask should cover
%               either the heart, muscular, or surrounding tissue to weight
%               the lung image

%               lung_val: this is the assigned mask label from
%               Amira/ITK_Snap segmentations. This mask should cover the
%               entire lung vol (assuming image reconstructed at
%               expiration)

%               savefile: string that should be the filename you save the
%               output structure variable, ImgResults

% Output:       ImgResults: Structure variable containing each metric
%               computated

% Additional code required/suggested: DICOMLoad.m      | TIFFLoad.m | imslice.m
%                                     calcImageSlice.m | Volume.m   | showImageSlice.m
%                                     LoadData_Amira.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Changing directory to working folder containing images and masks to load in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
cd(path);

% Loading in End-Expiration DICOM binary masks
disp('Load all Exp DICOM files containing masks generated from Amira: ');
exp_masks = DICOMLoad;

% Loading in End-Inspiration DICOM binary masks
% Note: inspiratory masks here should only have a single label
disp('Load all Insp DICOM files containing masks generated from Amira: ');
insp_masks = DICOMLoad;
insp_masks(insp_masks > 0) = 1;

% Loading in TIFF image slices
disp('Load all TIFF files containing MRI images of mouse: ');
mri_img = TIFFLoad;

% % Alternative way to load in both masks and images
% % 1. By LoadData_Amira.m
% [file, base] = uigetfile('*.am', 'Select the segmented mask file generated in Amira');
% amira_masks = double(LoadData_Amira(fullfile(base,file)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generation of binary masks from the loaded in DICOM files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating binary lung mask from all amira masks
lung_mask = exp_masks;
lung_mask(lung_mask ~= lung_val) = 0;
lung_mask(lung_mask > 0) = 1;   

% Creating binary heart (Weighting ROI) mask from all amira masks
roi_mask = exp_masks;
roi_mask(roi_mask ~= roi_val) = 0;
roi_mask(roi_mask > 0) = 1;

% Making the lung and weighting roi images
lung_img  = lung_mask .* mri_img;
roi_img = roi_mask .* mri_img;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the mean of the weighting ROI image (i.e., heart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
roi_mean = mean(roi_img(roi_img > 0));

% Weighting lung image with ROI mean
wtLung_img = lung_img / roi_mean;

% Display the weighted lung image 
imslice(wtLung_img, 'Weighted lung img from ROI mean');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Using baseline threshold from parenchyma segmented images at Day 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contact corresponding author if you wish to see this data that was used 
% to generate this threshold

base_mean = 0.5899;
base_std  = 0.0727;
factor    = 3;

base_thres = base_mean + (factor * base_std);

% Creating the thresholded image vector
wtLung_vector = wtLung_img(wtLung_img > base_thres);

% Creating the thresholded, and weighted, image containing only high signal 
highSig_img = wtLung_img;
highSig_img(highSig_img < base_thres) = 0;

% Display the thresholded image of weighted high signal volume
figure; imslice(highSig_img, 'Thresholded lung img after weighting');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating metrics from mouse images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean and coefficient of variation of the weighted lung image
wtLung_mean = mean(wtLung_img(lung_mask == 1));
wtLung_std  = std(wtLung_img(lung_mask == 1));
wtLung_cv   = wtLung_std / wtLung_mean;

% Percentage of high signal volume (PHSV)
% Note: 3.125e-5 is calculated from FOV / Matrix Size
% e.g, FOV = 3.2 x 3.2 x 6.4 (cm^3) / Matrix Size = (128 x 128 x 128) = 3.125e-5
voxvol = (3.2*3.2*6.4) / (128^3);
lungvol_exp = sum(lung_mask(lung_mask ==  1)) * voxvol; % lung volume in mL
highsig_vol = length(wtLung_vector) * voxvol;
PHSV        = (highsig_vol / lungvol_exp) * 100; % expressed as a percentage

% Percent tidal volume (PTV), weighted by inspiratory volume
lungvol_insp = sum(insp_masks(insp_masks == 1)) * voxvol;

tidalvol = (lungvol_insp - lungvol_exp); % expressed in mL
PTV = ((lungvol_insp - lungvol_exp) / lungvol_insp) * 100;% expressed as a percentage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saving the output variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImgResults.LungMean = wtLung_mean;
ImgResults.LungCV   = wtLung_cv;
ImgResults.PHSV     = PHSV;
ImgResults.PTV      = PTV;
ImgResults.ROIMean  = roi_mean;
ImgResults.LungVol.ExpVol   = lungvol_exp;
ImgResults.LungVol.InspVol  = lungvol_insp;
ImgResults.LungVol.TidalVol = tidalvol;
ImgResults.Threshold.BaseMean = base_mean;
ImgResults.Threshold.BaseSTD  = base_std;
ImgResults.Threshold.Factor   = factor;
ImgResults.Images.RawImg = mri_img;
ImgResults.Images.ThresImg = highSig_img;
ImgResults.Masks.AmiraMasks = exp_masks;
ImgResults.Masks.LungExpVol = lung_mask;
ImgResults.Masks.LungInspVol = insp_masks;
save(savefile, 'ImgResults');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


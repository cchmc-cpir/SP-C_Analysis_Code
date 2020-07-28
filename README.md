# SP-C_Analysis_Code
Image reconstruction and analysis code used for repeatd mesures on SP-C mice mutants

# Function Descriptions
Bruker_Load.m: Loads in raw free induction decay (FID) file into Matlab workspace as complex double

DICOMLoad.m: Loads in DICOM (.dcm) files of manually segmeneted, binary masks exported from Amira software. The user selects the folder containing the .dcm files and all files are automatically selected.

ImgAnalysis_SPC.m: Analysis function that measures MR-derived metrics (mean lung signal, lung coefficient of variation, tidal volume, and percentage of high signal volume). There are 4 input arguments that must be fullfilled:
  1. path: the directory that houses all image data and binary masks
  2. roi_val: label value assigned from segmentation software (e.g., Amira, ITK-Snap) pertaining to a soft tissue region (i.e., heart) to weight the segmented lung image.
  3. lung_val: label value assigned from segmetnation software pertaining to the lung volume.
  4. savefile: string variable name that saves output results as a Matlab data file.
  
  LoadData_Amira.m: Alternative means to load in segmented binary masks. If you use Amira, this function will read in the Amira (.am) file is the fullfile path is specified.
  
  TiffLoad.m: Loads in TIFF (.tif) files of reconstructed image that was retrospectively gated at end-expiration. The user selects the folder containing the .tif files and all files are automatically selected.
  
quick_retrogating_file.m: Performs retrospective gating on a loaded in FID file (must be complex variable) that is reshaped in the proper format (e.g., 81 x 44,334; where 81 is the number of FID points and 44,334 is the number of radial projections). The FID is gated with default values set for end-expriation/inspiration, however, the user determines if these parameters are sufficient or can select and choose new gating parameters.
  
read_UTEmethod.m: Reads in and extracts scan parameters like number of FID points, radial projections, FID, etc. The folder path containing the method file must be specified to run this function.

% could wait for 2nd round to convert to python, before that could use hard
% code in python
% [iField,voxel_size,matrix_size,CF,delta

%[TE,TE,B0_dir] = Read_DICOM(get(handles.DataLocationEditable,'string'));
load('medi_siemens_data.mat')

iMag = sqrt(sum(abs(iField).^2,4));
% need convert to python
[iFreq_raw N_std] = Fit_ppm_complex(iField);
save('matlabdata_iFreq_raw.mat','iFreq_raw')
% could wait for 2nd round to convert to python
% region growing method
% iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);

% Spatial phase unwrapping %%%
% if large fringe lines persists, try 
% need  convert to python
iFreq = unwrapLaplacian(iFreq_raw, matrix_size, voxel_size);
save('matlabdata_iFreq.mat','iFreq')
% Mask = genMask(iField, voxel_size);
% use exe from FSL instead, no need to convert to python
Mask = BET(iMag, matrix_size, voxel_size);
% qsm_filepath = 'D:\eclipse-workspace\OncoImageAnalysis\src\QSM\data\004_QSM_MONO_8TE_IPAT2_68phpf_NoFC_mask.nii';
% qsm_filepath_mask = 'D:\eclipse-workspace\OncoImageAnalysis\src\QSM\data\004_QSM_MONO_8TE_IPAT2_68phpf_NoFC_matmask.nii';
% yx_save_nii_file(Mask(:,end:-1:1,:),qsm_filepath,qsm_filepath_mask)

% need convert to python
R2s = arlo(TE, abs(iField));

% need convert to python


Mask_CSF = extract_CSF(R2s, Mask, voxel_size);

save('matlabdata_iFreq_FOR RDF.mat','iFreq')
save('matlabdata_Mask_FOR RDF.mat','Mask')
save('matlabdata_N_std_FOR RDF.mat','N_std')
save('matlabdata_matrix_size_FOR RDF.mat','matrix_size')
save('matlabdata_voxel_size_FOR RDF.mat','voxel_size')
save('matlabdata_ B0_dir_FOR RDF.mat', 'B0_dir')
% need convert to python
RDF = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir);


figure,imshow(RDF(:,:,30),[-1,1]);
%%%% Background field removal using Projection onto Dipole Fields
%%%% NMR Biomed 2011;24(9):1129-36.
%%%% MRM 2010;63(1):194-206

% Zhou et al. NMR in Biomed 27 (3), 312-319, 2014
% could wait for 2nd round to convert to python
% RDF = LBV(iFreq,Mask,matrix_size,voxel_size);

% Before running MEDI, variables need to be saved
%save RDF.mat RDF iFreq iFreq_raw iMag N_std Mask matrix_size...
    % voxel_size delta_TE CF B0_dir Mask_CSF;
% need to convert to python, slow speed, better use cython
QSM = MEDI_L1('lambda',1000,'percentage',0.9, 'smv',5);
figure,imshow(QSM(:,:,30),[])
function scanner_coord_raw_native = get_scanner_coords_longitudinal(S1_filename, coords, S2_filename)
% usage: scanner_coords_input_signed = get_scanner_coords_jong('<nifti from S1>', <xyz coords from S1>, '<nifti from S2>');
%
% Inputs:
% S1_filename = t1w image from session S1
% S2_filename = t1w image from session S2
% coords = in S1 physical space,  is the point that we want to transform
% to the S2 brain in scanner space in the following order: [R/L, A/P, S/I]
% e.g., if on the scanner the voxel center in L-DLPFC is [S30, L22, A28],
% then the input coords = [-22,28,30]
%
% Outputs:
% scanner_coords_raw_native = coordinates to be entered into the scanner

addpath(genpath('/home/cniuser/git/spm12/'));
rmpath(genpath('/home/cniuser/git/spm12/external/fieldtrip/'));

if ~exist('S2_filename', 'var') || isempty(S2_filename), S2_filename = 'anat.nii.gz'; end

S1 = spm_vol(S1_filename);
S2 = spm_vol(S2_filename);

spm_figure('CreateWin', 'Graphics');
x = spm_coreg(S1,S2)

M = spm_matrix(x);

coords_to_be_transformed = [coords 1]';
fourd_coord = mtimes(M,coords_to_be_transformed)';
scanner_coord_raw_native = fourd_coord(1:3);

% convert coordinates to scanner input 
% ?? i.e., flip around, +Z to S (first), +X to R (second), +Y to A (third)
scanner_coord_input_reordered = NaN(1:3);
scanner_coord_input_reordered(1) = scanner_coord_raw_native(3);
scanner_coord_input_reordered(2) = scanner_coord_raw_native(1);
scanner_coord_input_reordered(3) = scanner_coord_raw_native(2);
scanner_coord_input_reordered_abs = abs(scanner_coord_input_reordered);

coord_signs = sign(scanner_coord_input_reordered);

if coord_signs(1) == 1
    coord_1_orient = 'S';
else
    coord_1_orient = 'I';
end

if coord_signs(2) == 1
    coord_2_orient = 'R';
else
    coord_2_orient = 'L';
end

if coord_signs(3) == 1
    coord_3_orient = 'A';
else
    coord_3_orient = 'P';
end

fprintf(['scanner coordinate: \n  ' coord_1_orient num2str(roundn(scanner_coord_input_reordered_abs(1),-1)) ' ' coord_2_orient num2str(roundn(scanner_coord_input_reordered_abs(2),-1)) ' ' coord_3_orient num2str(roundn(scanner_coord_input_reordered_abs(3),-1)) ' \n\n'])








function scanner_coords_raw_native = get_scanner_coords_MNI(labels, coords, individual_nii, template, individual_anat, check_xf)
% usage: scanner_coords_input_signed = get_scanner_coords_MNI(labels, coords, individual_nii, template, check_xf);
%
% NOTE: mni_template_coords = in physical space, and are the points that we
% want to be transformed to the individuals brain in scanner space, should
% be a matrix with three columns and a row for each coordinate
% coordinates from AFNI need to have the first two coordinates flipped
% (from the coordinates in the upper left of the viewer box under [order:
% RAI=DICOM]) to match the physical space coordinates in the second column
% of the FSL viewer coordinates, which is what this program takes
% approximate coords and boxes:
%
% the output coordinates: orient = {'R','A','S';'L','P','I'}'; - - - + + +
% fsl: RL = +-; AP = +-; SI = +-
%
% Inputs:
%   labels              label for the voxels of interest. 1xn cell
%   coords              voxel coordinates in the MNI space. nx3 array
%   individual_nii      individual's T1w image. Default is anat.nii.gz
%   template            template nifti in which the coordinates were drawn, e.g. masked MNI brain.
%                       Default is MNI152_T1_2mm_brain.nii
%   check_xf            flag for showing SPM normalization results. Default is true
%
% Outputs:
%   scanner_coords      coordinates to be entered into the scanner


% setup environment and variables
addpath(genpath('~/git/spm12'));
setenv('FSLDIR', '/etc/fsl/5.0');
fsldir = getenv('FSLDIR');
if ~exist('individual_nii', 'var') || isempty(individual_nii), individual_nii = 'anat.nii.gz'; end
if ~exist('template', 'var') || isempty(template), template = 'MNI152_T1_2mm_brain.nii'; end
if ~exist('check_xf', 'var') || isempty(check_xf), check_xf = true; end

% default ROI in MNI space
if ~exist('labels', 'var') || isempty(labels), labels = {'lDLPFC'}; end
if ~exist('coords', 'var') || isempty(coords), coords = [-29,26,32]; end

mni_template_coords = coords';

% mask individual T1 and uncompress nifti because spm requires .nii
individual_nii_masked = [individual_nii(1:end-7) '_brain.nii.gz'];
cmd = ['bet ' individual_nii ' ' individual_nii_masked ' -f 0.4'];
command = sprintf('/bin/bash -c ''. %s/fsl.sh; fsl5.0-%s''', fsldir, cmd);
status = system(command);
individual_nii_masked_unzip = gunzip(individual_nii_masked);
individual_nii_masked_unzip = individual_nii_masked_unzip{1};

% calculate coordinates
ni = spm_vol_nifti(individual_nii_masked_unzip);
spm_defaults; global defaults;
params = defaults.old.normalise.estimate; 
if check_xf, params.graphics = 1; spm_figure('Create','Graphics'); end  % if true, show the normalization figures
sn = spm_normalise(template, ni, '', '', '', params);
scanner_coords_raw_native = getScannerCoordsFromSn(sn, mni_template_coords)';

% convert coordinates to scanner input 
% i.e., flip around, +Z to S (first), +X to R (second), +Y to A (third)
scanner_coords_input_reordered = NaN(size(scanner_coords_raw_native));
scanner_coords_input_reordered(:,1) = scanner_coords_raw_native(:,3);
scanner_coords_input_reordered(:,2) = scanner_coords_raw_native(:,1);
scanner_coords_input_reordered(:,3) = scanner_coords_raw_native(:,2);
scanner_coords_input_reordered_abs = abs(scanner_coords_input_reordered);

coord_signs = sign(scanner_coords_input_reordered);

fprintf('\nScanner input coords: \n\n')

for region_num = 1:length(labels)

    if coord_signs(region_num,1) == 1
        coord_1_orient = 'S';
    else
        coord_1_orient = 'I';
    end
    
    if coord_signs(region_num,2) == 1
        coord_2_orient = 'R';
    else
        coord_2_orient = 'L';
    end
    
    if coord_signs(region_num,3) == 1
        coord_3_orient = 'A';
    else
        coord_3_orient = 'P';
    end
    
    fprintf([cell2mat(labels(region_num)) ': \n   ' ...
        coord_1_orient num2str(roundn(scanner_coords_input_reordered_abs(region_num,1),-1)) ' ' ...
        coord_2_orient num2str(roundn(scanner_coords_input_reordered_abs(region_num,2),-1)) ' ' ...
        coord_3_orient num2str(roundn(scanner_coords_input_reordered_abs(region_num,3),-1)) ' \n\n'])

end

 % check transform
%  outIm = mrAnatResliceSpm(ni.img, sn, [-80,-120,-60;80,120,120]);
%  showMontage(outIm);

delete(individual_nii_masked, individual_nii_masked_unzip);  % delete intermediate files

return;


%_________________________________________________________________________________
% code from mrAnatGetImageCoordsFromSn in vistasoft
function scanner_coords = getScannerCoordsFromSn(sn, template_coords)
    
[x,y,z] = deal(template_coords(1,:),template_coords(2,:),template_coords(3,:));
% from template physical space to voxel space
Mult = inv(sn.VG.mat);  
[x,y,z] = mmult(x,y,z,Mult);
% from template to scanner space 
Mult = sn.VF.mat*sn.Affine;  % affine 
if ~isempty(sn.Tr)    % add nonlinear components 
    dim = sn.VG.dim;
    if(length(dim)==3), dim(4) = 0; end
    [x,y,z] = build_transform2(sn.Tr, [size(sn.Tr); dim], x, y, z);
end
[x,y,z]  = mmult(x, y, z, Mult);

scanner_coords = [x;y;z];

return;

function [TX,TY,TZ] = build_transform2(T,dim,TX,TY,TZ)
BX = basis_funk(TX,dim(2,1),dim(1,1)); 
BY = basis_funk(TY,dim(2,2),dim(1,2));
BZ = basis_funk(TZ,dim(2,3),dim(1,3));
for i3=1:dim(1,3)
    for i2=1:dim(1,2)
		B2 = BZ(:,:,i3).*BY(:,:,i2);
        for i1=1:dim(1,1)
			B  = B2.*BX(:,:,i1);
			TX = TX + T(i1,i2,i3,1)*B;
			TY = TY + T(i1,i2,i3,2)*B;
			TZ = TZ + T(i1,i2,i3,3)*B;
        end
    end
end
return;

function B = basis_funk(X,N,kk)
B = zeros([size(X) kk]);
B(:,:,1) = ones(size(X))/sqrt(N);
for k=2:kk
    B(:,:,k) = sqrt(2/N)*cos((X-0.5)*(pi*(k-1)/N));
end
return;

function [X2,Y2,Z2] = mmult(X1,Y1,Z1,Mult)
if length(Z1) == 1
	X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + (Mult(1,3)*Z1 + Mult(1,4));
	Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + (Mult(2,3)*Z1 + Mult(2,4));
	Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + (Mult(3,3)*Z1 + Mult(3,4));
else
	X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
	Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
	Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
end
return;

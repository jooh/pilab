
% replace this value with the directory that meets the above requirements
resdir = pwd;
% construct the Volume instances in one line. epivol contains in-mask data, ignoring all
% voxels with NaNs at any time point (these are usually introduced by motion correction).
[epivol,designvol] = spm2vol(fullfile(resdir,'SPM.mat'),'mask',fullfile(resdir,'pilab_mask.nii'));
% it's worth taking a closer look at epivol, which contains a few bells and whistles we
% didn't see in tutorial 1
epivol

% why the vector of ones? Well, we just want to write out a volume of where every in-mask voxel is 1
% (and everything else gets left to the default value of zero)
data2file(epivol,ones(1,epivol.nfeatures),fullfile(resdir,'pilab_mask_updated.nii'));

% here the input is the first row of the data matrix, that is, all the voxels from the first TR.
% makeimagestack is a nice function from Kendrick Kay for tiling slices through a 3D volume to a
% 2D plane for visualisation.
imagesc(makeimagestack(data2mat(epivol,epivol.data(1,:))));
% pedantic styling below
axis('off');
set(gca,'dataaspectratio',[epivol.voxsize(1:2) 1])
chandle = colorbar;
set(chandle,'position',get(chandle,'position') .* [1.1 1 .5 .5])

% add the regions of interest
roivol = roidir2vol(fullfile(resdir,'ROIs'))

[rfepi,rfroi] = intersectvols(epivol,roivol);
rfepi

roiepi = rfepi(:,rfroi.data(1,:)~=0)

imagesc(makeimagestack(data2mat(roiepi,roiepi.data(1,:)))~=0);
axis('off');
set(gca,'dataaspectratio',[epivol.voxsize(1:2) 1])

% heavy lifting ahead - check for matlabpool availability and start a pool if one isn't running already
if ~isempty(which('matlabpool')) && ~matlabpool('size')
    matlabpool local
end
% passing the epivol instance means we run the searchlight mapping on the same voxels as the main dataset.
% You could also enter a nifti instead, e.g. the file we generated above (pilab_mask_updated.nii)
[searchrois,diagnostic] = mask2searchrois(epivol,'nvox',100);

% visualise the radius of each searchlight
imagesc(makeimagestack(data2mat(epivol,diagnostic.r)));
axis('off');
set(gca,'dataaspectratio',[epivol.voxsize(1:2) 1])
chandle = colorbar;
title('radius (mm) per searchlight')
% typical radius
mode(diagnostic.r)

% visualise how often a each voxel is sampled by a searchlight
imagesc(makeimagestack(data2mat(epivol,diagnostic.nsampled)));
axis('off');
set(gca,'dataaspectratio',[epivol.voxsize(1:2) 1])
chandle = colorbar;
title('number of searchlights per voxel')
% typical number of searchlights per voxel
mode(diagnostic.nsampled)

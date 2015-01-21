% wrapper for Searchlight. Calculate all in-mask searchlight for some
% maskpath and return an MriVolume instance of the mapped searchlight rois.
% 
% INPUTS
% maskpath: brain mask (path to nifti or MriVolume instance)
% mapmode: radius or nvox
% radvox: parameter for mapmode (mm radius or number of voxels)
%
% OUTPUTS
% roivol: MriVolume instance with one row per searchlight and one column
%   per in-mask feature. NB, roivol.data will be sparse to conserve memory.
% diagnostic: struct with information about each searchlight: its radius
%   (r), number of voxels (nsphere), and number of searchlights that
%   samples that voxel (nsampled).
%
% NB! diagnostic.nsphere and diagnostic.nsampled will be identical for
% radius-based mapping.  This is because if a voxel x1 is outside a
% searchlight centered on x2, then a searchlight on x1 necessarily doesn't
% include x1, and conversely. So the distribution of nvoxels == the
% distribution of nsamples or intuitively, a searchlight that includes many
% voxels will also be included as a voxel in many other searchlights. But
% this does not hold for nvox-based mapping.
%
% [roivol,diagnostic] = mask2searchrois(maskpath,mapmode,radvox)
function [roivol,diagnostic] = mask2searchrois(maskpath,mapmode,radvox)

sl = Searchlight(maskpath,mapmode,radvox);

% new sparse formulation - store radius information directly in roivol
% instance.
spheres = sparse(sl.vol.nfeatures,sl.vol.nfeatures);
r = NaN([1 sl.vol.nfeatures]);

% run
fprintf('mapping %d searchlights...\n',sl.vol.nfeatures);
tic;
parfor n = 1:sl.vol.nfeatures
    % get sphere index
    % (some extra dribbling is necessary here to avoid confusing
    % parfor)
    sp = sparse(1,sl.vol.nfeatures);
    inds = sl.mapinds(n);
    sp(inds) = sl.distances(inds);
    spheres(n,:) = sp;
    r(n) = sl.radius;
end
fprintf('finished in %s.\n',seconds2str(toc));
diagnostic.r = r;
% number of voxels in each sphere
sb = spheres;
% make sure sum operation below counts each sphere/voxel once
sb(sb~=0) = 1;
diagnostic.nsphere = full(sum(sb,2)');
% number of spheres that sampled each voxel
diagnostic.nsampled = full(sum(sb,1));
% preserve class
roivol = feval(class(sl.vol),spheres,sl.vol);

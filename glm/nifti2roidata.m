% convert a set of nifti volumes to the 'roidata' structure used elsewhere.
%
% INPUTS:
% niftis: char array of niftis (one per contrast)
% mask: logical 3D of (group) in-mask voxels
% targetfield: what to name the roidata field that contains the actual data
% (default 't')
%
% roidata = nifti2roidata(niftis,mask,targetfield)
function roidata = nifti2roidata(niftis,mask,targetfield,setclass)

if ieNotDefined('targetfield')
    targetfield = 't';
end

if ieNotDefined('setclass')
    setclass = 'double';
end

V = spm_vol(niftis);
ncon = numel(V);

% number of valid entries
nroi = sum(mask(:));

roidata.(targetfield) = NaN([ncon nroi],setclass);
roidata.cols_roi = cell(1,nroi);
roidata.rows_contrast = cell(ncon,1);

for n = 1:ncon
    [~,roidata.rows_contrast{n},~] = fileparts(V(n).fname);
    xyz = spm_read_vols(V(n));
    roidata.(targetfield)(n,:) = xyz(mask);
end

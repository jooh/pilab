% variance smooth permutation test using 'pseudo' T statistics (see e.g.
% Nichols & Holmes, 2001).
% 
% INPUTS:
% permest: permutation distribution of effect (2d with permutation in
%   columns, true perm is first column).
% permvar: permutation distribution of variance estimate
% mask: path to SPM nifti or 3D logical matrix (if mat, must define voxsize
%   named input below)
%
% NAMED INPUTS:
% fwhm: (8) in mm (or whatever unit voxsize is in)
% voxsize: ([])
% tail: (default 'right')
%
% [p,pseudot] = permpfwe_varsmooth(permest,permvar,mask,varargin)
function [p,pseudot] = permpfwe_varsmooth(permest,permvar,mask,varargin)

getArgs(varargin,{'fwhm',10,'voxsize',[],'tail','right'});
assert(ndims(permest)==2 && ndims(permvar)==2,...
    'inputs must be 2D with permutations in columns');

% parse mask input (nb ismatrix is no good here because this returns false
% for logical)
if ismat(mask)
    assert(~isempty(voxsize),['voxsize must be defined if entering ' ...
        'mask as matrix rather than path to nifti']);
else
    V = spm_vol(mask);
    mask = spm_read_vols(V) > 0;
    voxsize = vox2mm(V);
end

% smooth the variance estimate
smoothvar = smoothdatavecs(permvar',fwhm,mask,voxsize);

% get the pseudo T - null effect estimate / smoothed variance estimate
pseudotdist = permest' ./ smoothvar;

% max stat p value as usual
p = permpfwe(pseudotdist,tail);
pseudot = pseudotdist(1,:);

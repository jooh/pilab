% Convenience function for loading in-mask voxels for a set of volumes
% using spm_get_data.
%
% Inputs:
% paths: cell or char array of file paths to SPM-readable volumes
% mask: 3D matrix in same shape as volumes
%
% If you specify the mask output this function will strip out voxels that
% are 'bad' (ie, NaN or 0 at any time point) and return both cleaned up
% data and a cleaned up mask. If you specify only data output we assume you
% don't want any cleaning of data or mask.
%
% [data,[mask]] = loadmaskedvolumes(paths,mask)
function [data,mask] = loadmaskedvolumes(paths,mask)

% ensure logical
mask = mask > 0;
maskind = find(mask);
[x,y,z] = ind2sub(size(mask),maskind);

% 2d data (time by voxels)
% note that the final false flag disables internal SPM error checking.
% This speeds things up but does assume that your mask is sensible.
data = spm_get_data(spm_vol(paths),[x y z]',false);

% if you asked for the mask back, assume you want to strip bad data too
if nargout > 1
    badinds = any(isnan(data) | data==0,1);
    data(:,badinds) = [];
    mask(maskind(badinds)) = false;
    if any(badinds)
        logstr('removed %d voxels (%2.0f%% of total)\n',...
            sum(badinds),100*sum(badinds) / numel(badinds));
    end
end

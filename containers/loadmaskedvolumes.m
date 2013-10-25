% Convenience function for loading in-mask voxels for a set of volumes.
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
data = spm_get_data(spm_vol(paths),[x y z]');

% if you asked for the mask back, assume you want to strip bad data too
if nargout > 1
    badinds = any(isnan(data) | data==0,2);
    data(:,badinds) = [];
    mask(maskind(badinds)) = false;
    if any(badinds)
        fprintf('removed %d voxels (%d percent of total)\n',...
            sum(badinds),100 * sum(badinds) / length(badinds));
    end
end

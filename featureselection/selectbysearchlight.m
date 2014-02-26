% return logical feature indices based on selecting searchlights according
% to some (independent) result.
%
% INPUTS (roivol and result are assumed to be in register)
% roivol: searchlight volume most likely. 
% result: some effect estimate in vector, matrix or path to nifti form
% n: number of features to use for selection
%
% OPTIONAL, NAMED INPUTS (default):
% sortdirection: (ascend)
% abssort: (false) abs transform the result before ranking
% selectmode: (center) return the winning n centers, alternatively
%   union - any feature in one of the selected searchlights
%   intersect - any feature in all the selected searchlights
%
% OUTPUTS:
% inds - logical or numeric indices
% thresh - result statistic at the threshold
%
% [inds,thresh] = selectbysearchlight(roivol,result,n,[varargin])
function [inds,thresh] = selectbysearchlight(rois,result,n,varargin)

getArgs(varargin,{'sortdirection','descend','abssort',false,...
    'selectmode','centre'});

% parse the result to obtain a column vector for ranking
if ischar(result)
    result = spm_read_vols(spm_vol(result));
end
switch ndims(result)
    case 2
        assert(isrow(result) || iscolumn(result),'2d input must be vector')
        % as column vector
        result = result(:)';
    case 3
        result = result(rois.mask);
    otherwise
        error('could not parse result dimensionality: %d',ndims(result));
end
if abssort
    result = abs(result);
end

% rank and obtain indices for feature selection
[resv,resind] = sort(result,sortdirection);
switch selectmode
    case 'center'
        % nice and easy
        inds = ind2logical(size(result),resind(1:n));
        thresh = resv(n);
    case {'union','intersect'}
        if strcmp(selectmode,'union')
            indfun = @any;
        else
            indfun = @all;
        end
        nspheres = 0;
        inds = [];
        % keep adding searchlights until we reach the desired total ROI
        % size (err on the size of over-shooting)
        while sum(inds)<n
            nspheres = nspheres + 1;
            inds = indfun(full(rois.data(resind(1:nspheres),:))~=0,1);
            thresh = resv(nspheres);
        end
    otherwise
        error('unknown selectmode: %s',selectmode)
end

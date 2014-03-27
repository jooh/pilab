% average a set of identically-shaped matrices.
%
% result = matmean(varargin);
function result = matmean(varargin);

if ~nargin
    result = [];
    return
end

catdim = ndims(varargin{1})+1;
result = mean(cat(catdim,varargin{:}),catdim);

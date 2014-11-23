% average a set of identically-shaped matrices.
%
% result = matmean(varargin);
function result = matmean(varargin);

switch nargin
    case 0
        result = [];
    case 1
        result = varargin{1};
    otherwise
        catdim = ndims(varargin{1})+1;
        result = mean(cat(catdim,varargin{:}),catdim);
end

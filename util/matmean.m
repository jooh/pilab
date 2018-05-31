% average a set of identically-shaped matrices.
%
% result = matmean(varargin);
function result = matmean(varargin);

result = matfun(@mean,varargin{:});

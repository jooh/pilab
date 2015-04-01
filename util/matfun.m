% apply some function to the identically-shaped input varargin. The funhand
% must take an NDarray as its first input and a dim argument as its second
% (e.g., mean, median). To work with other functions, try anonymous
% functions, e.g.  @(x,dim)std(x,[],dim) to get around the fact that std
% uses the third argument for dim.
%
% result = matfun(funhand,varargin);
function result = matfun(funhand,varargin);

switch nargin
    case 0
        result = [];
    otherwise
        catdim = ndims(varargin{1})+1;
        result = feval(funhand,cat(catdim,varargin{:}),catdim);
end

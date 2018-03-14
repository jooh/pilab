% You know how Matlab doesn't have a built-in function for calculating the
% standard error of the sample mean?
%
% y = sterr(data,dim);
function y = sterr(data,dim);

if ~exist('dim','var') || isempty(dim)
    dim = 1;
end

y = std(data,[],dim) ./ sqrt(size(data,dim));

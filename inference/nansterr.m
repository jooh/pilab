%
% y = nansterr(data,dim);
function y = nansterr(data,dim);

if ieNotDefined('dim')
    dim = 1;
end

y = nanstd(data,[],dim) ./ sqrt(sum(~isnan(data),dim));

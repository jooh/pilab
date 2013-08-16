% Conduct a percentile bootstrap on the distribution in bootest (with
% boostraps in the final dimension). return the median estimates and
% standard errors (half difference between 16th and 84th percentiles).
%
% [estimates,sterrs] = bootprctile(bootest)
function [estimates,sterrs] = bootprctile(bootest)

% 2nd or 3rd I would expect
bootdim = ndims(bootest);
percs = prctile(bootest,[16 50 84],bootdim);
switch bootdim
    case 2
        estimates = percs(:,2);
        sterrs = diff(percs(:,[1 3]),1,2)/2;
    case 3
        estimates = percs(:,:,2);
        sterrs = diff(percs(:,:,[1 3]),1,3)/2;
    otherwise
        error('no bootprctile support for %d inputs',bootdim)
end

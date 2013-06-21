% Use a bootstrap distribution ( by nboot)
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
        error('bootprctile not support for %d inputs',bootdim)
end

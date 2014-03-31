% return the p values for a permutation distribution (with permutations in
% the final dimension). We assume that the first entry is the true sample.
% tail is 'right' (default), 'left' or 'both' for 2-tailed inference.
%
% p = permpvalue(permdist,tail)
function p = permpvalue(permdist,tail)

if ieNotDefined('tail')
    tail = 'right';
end

permdim = ndims(permdist);
datadims = size(permdist);

% there must be a more general way
switch permdim
    case 2
        truerep = repmat(permdist(:,1),[1 datadims(permdim)]);
    case 3
        truerep = repmat(permdist(:,:,1),[1 1 datadims(permdim)]);
    case 4
        truerep = repmat(permdist(:,:,:,1),[1 1 1 datadims(permdim)]);
    otherwise
        error('no permpvalue support for %d inputs',permdim);
end

switch tail
    case 'right'
        ngreater = sum(permdist >= truerep,permdim);
    case 'left'
        ngreater = sum(permdist <= truerep,permdim);
    case 'both'
        % two-tailed - absolute value (ie pos or neg) exceeding the
        % absolute nulls
        ngreater = sum(abs(permdist) >= abs(truerep),permdim);
    otherwise
        error('unknown tail: %s',tail);
end

% adjust for nan perms
nvalid = datadims(permdim) - sum(isnan(permdist),permdim);
p = ngreater ./ nvalid;
% if nvalid==0 you get inf above
p(isinf(p)) = NaN;

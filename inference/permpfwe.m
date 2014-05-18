% Return the nichols / holmes p value for a null distribution
% permutations are in rows, features in columns. The first row is assumed
% to contain the true, unpermuted effect.
%
% [pfwe,pthresh] = maxstatpfwe(nulldist,tail)
function [pfwe,pthresh] = maxstatpfwe(nulldist,tail)

[nperms,ndata] = size(nulldist);
assert(ndims(nulldist)==2,'can only support 2d inputs, got %dd',...
    ndims(nulldist));

if ieNotDefined('tail')
    tail = 'right';
end

switch tail
    case 'right'
        % do nothing
    case 'left'
        nulldist = nulldist * -1;
    case 'both'
        nulldist = abs(nulldist);
    otherwise
        error('unknown tail: %s',tail);
end

% assume truestat is first entry in nulldist
truestat = nulldist(1,:);

% distribution of max statistics across permutations
maxstats = max(nulldist,[],2);

% the FWE p value is the percentile of truestat in the max distribution
% (note btw that the uncorrected p value never enters into this
% calculation)
pfwe = sum(bsxfun(@ge,maxstats,truestat),1) / nperms;

pthresh = prctile(maxstats,95);

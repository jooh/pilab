% Return the nichols / holmes p value for a null distribution.  features
% are in rows, permutations are in columns. The first column is assumed to
% contain the true, unpermuted effect.
%
% [pfwe,pthresh] = permpfwe(nulldist,tail)
function [pfwe,pthresh] = permpfwe(nulldist,tail)

[ndata,nperms] = size(nulldist);
assert(ismatrix(nulldist),'can only support 2d inputs, got %dd',...
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
    case 'none'
        pfwe = NaN([ndata 1]);
        pthresh = NaN([ndata 1]);
        return;
    otherwise
        error('unknown tail: %s',tail);
end

% assume truestat is first entry in nulldist
truestat = nulldist(:,1);

% distribution of max statistics across permutations
maxstats = max(nulldist,[],1);

% the FWE p value is the percentile of truestat in the max distribution
% (note btw that the uncorrected p value never enters into this
% calculation)
pfwe = sum(bsxfun(@ge,maxstats,truestat),2) / nperms;

pthresh = prctile(maxstats,95);

% Return the nichols / holmes p value for a null distribution
% pfwe = maxstatpfwe(nulldist)
function pfwe = maxstatpfwe(nulldist)

[nperms,ndata] = size(nulldist);

% assume truestat is first entry in nulldist
truestat = nulldist(1,:);

% distribution of max statistics across permutations
maxstats = max(nulldist,[],2);

% the FWE p value is the percentile of truestat in the max distribution
% (note btw that the uncorrected p value never enters into this
% calculation)
pfwe = sum(repmat(maxstats,[1 ndata]) > repmat(truestat,[nperms 1]),1) ...
    / nperms;

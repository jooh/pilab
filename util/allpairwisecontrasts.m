% vectors for all pairwise contrasts given n conditions. For a version with
% labels, see roidata_allpairwisecontrasts.
%
% c = allpairwisecontrasts(n)
function c = allpairwisecontrasts(n)

% number of contrasts by number of regressors
outsize = [n nchoosek(n,2)];

% find indices for each pair
mask = tril(ones(n,n),-1);
[r,c] = find(mask);

% convert to linear indices into out matrix
pos = sub2ind(outsize,r,(1:outsize(2))');
neg = sub2ind(outsize,c,(1:outsize(2))');
c = zeros(outsize,class(n));
c(pos) = 1;
c(neg) = -1;
c = c';

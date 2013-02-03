% Vectorised form of pearson correlation coefficient for the special case
% where one input (column vector a) is to be correlated with a set of
% inputs (b), e.g. a predictor RDM correlated against each of a set of
% searchlight RDMs (each row is a different RDM).
%
% Uses \ on Z-scored data. Similar results can be had with corr(a,b) but
% this function is considerably faster.
% r = pearsonvec(a,b)
function r = pearsonvec(a,b)

[ar,ac] = size(a);
assert(ac==1,'a must be a column vector with one entry')
[br,bc] = size(b);
assert(ar==br,'a and b do not have the same number of rows')

% strip NaNs from a (ie, dissimilarities for which we have no prediction)
nans = isnan(a);
a(nans) = [];
b(nans,:) = [];

r = olsfit(zscore(a,0,1),zscore(b,0,1));

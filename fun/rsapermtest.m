% Permute rdm a and estimate a null distribution for its relationship with
% each of the rdms in b.
% nperms defaults to 1000
% distancemetric defaults to @spearmanvec
% [r,p,nulldist] = rsapermtest(a,b,[nperms],[distancemetric])
function [r,p,nulldist] = rsapermtest(a,b,[nperms],[distancemetric])

if ieNotDefined('nperms')
    nperms = 0;
end

if ieNotDefined('distancemetric')
    distancemetric = @spearmanvec;
end

a = asrdmmat(a);
b = asrdmmat(b);

% todo: get code for all permutations from old test. Revert to just
% randperming if the rdms are too big to avoid memory problems.
% also test that nperms is possible given rdm

% todo: parfor
for p = 1:nperms

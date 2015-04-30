% convert nchoosek(c,2) by n matrix to a c by c by n 3D stack of RDMs.
% Extends the builtin squareform by support for converting multiple vectors
% to RDMs in a single call.
%
% rdms = vec2rdm(vec)
function rdms = vec2rdm(vecs)

[npairs nrdm] = size(vecs);
assert(ismatrix(vecs),'input must be 2d with pairs in second dim')

% find the number of conditions given the number of pairs
n = npairs2n(npairs);
assert(n==round(n),'%d is not a valid number of pairs for RDM',npairs);

rdms = zeros([n,n,nrdm],class(vecs));
% fill in lower triangular vector
rdms(repmat(tril(true(n),-1),[1 1 nrdm])) = reshape(vecs,[1 npairs*nrdm]);
% and transpose to add lower triangular vector
rdms = rdms + permute(rdms,[2 1 3]);

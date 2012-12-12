% Convert a 3D stack of RDMs to vector form (lower triangular) with one
% row per RDM. This extends the Matlab built-in pdist by supporting
% vectorising multiple distance matrices (stacked in 3rd dim) in one
% operation.
% vecs = rdm2vec(rdms)
function vecs = rdm2vec(rdms)

[r,c,z] = size(rdms);

assert(r==c,'input must be square in first 2 dims!')

vecs = reshape(rdms(repmat(tril(true(r),-1),[1 1 z])),[nchoosek(r,2) z])';

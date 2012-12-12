% Convert a 3D stack of RDMs (or a struct array with an RDM field) to
% vector form (lower triangular) with one row per RDM. This extends the
% Matlab built-in pdist by supporting vectorising multiple distance
% matrices (stacked in 3rd dim) in one operation.
% vecs = rdm2vec(rdms)
function vecs = rdm2vec(rdms)

% unpack RSA toolbox-style rdm struct arrays
if isstruct(rdms)
    n = size(rdms(1).RDM,1);
    rdms = reshape([rdms.RDM],[n n length(rdms)]);
end

[r,c,z] = size(rdms);

assert(r==c,'input must be square in first 2 dims!')

vecs = reshape(rdms(repmat(tril(true(r),-1),[1 1 z])),[nchoosek(r,2) z])';

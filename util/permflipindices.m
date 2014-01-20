% generate nperms rows of 1:ncon logical indices, where 1 indicates a
% datapoint that should be flipped and 0 one that should be left alone.
% Used for a Nichols & Holmes-style sign flip permutation test.
%
% perminds = permflipindices(ncon,nperms)
function perminds = permflipindices(ncon,nperms)

npossible = 2^ncon;
if nperms > npossible
    warning(...
    'requested %d permutations but only %d are possible',nperms,npossible);
    nperms = npossible;
end

% this bit is a bit hairy (see Matlab's dec2bin for original
% implementation)
[~,e]=log2(npossible);
% perminds is now a logical matrix with perms in rows, and indices to flip
% in columns the true perm (all false) always ends up first
perminds = logical(rem(floor((npossible:-1:1)'*pow2(1-(e-1):0)),2));

% so this is exhaustive. If you wanted less than an exhaustive test...
if nperms < npossible
    % randomise order (add 1 to avoid original perm)
    randsamp = randperm(npossible-1) + 1;
    perminds = [perminds(1,:); perminds(randsamp(1:nperms-1),:)];
end

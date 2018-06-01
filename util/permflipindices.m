% generate nperms rows of 1:ncon logical indices, where 1 indicates a
% datapoint that should be flipped and 0 one that should be left alone.
% Used for a Nichols & Holmes-style sign flip permutation test.
%
% If ncon<26 we use an exhaustive method where we sample from all possible
% permutations. For larger ncon we use random number generation to approximate.
% Sampling is always done with the constraint that all permutations are unique.
%
% perminds = permflipindices(ncon,nperms)
function perminds = permflipindices(ncon,nperms)

% for reasons I don't entirely understand, you will only get the correct
% sampling proportions for the exhaustive version of you work in double rather
% than single precision.
ncon = double(ncon);
nperms = double(nperms);

npossible = 2^ncon;
if nperms > npossible
    warning(...
    'requested %d permutations but only %d are possible',nperms,npossible);
    nperms = npossible;
end

if ncon < 26
    % generate all possible permutations and then sample
    [~,e] = log2(npossible);
    perminds = logical(rem(floor((npossible:-1:1)'*pow2(1-(e-1):0)),2));
    if nperms < npossible
        % randomise order (add 1 to avoid original perm)
        randsamp = randperm(npossible-1) + 1;
        perminds = [perminds(1,:); perminds(randsamp(1:nperms-1),:)];
    end
else
    % need to sample from the random number generator to avoid running out of
    % memory (26 means we should stay under 4GB RAM, but careful here because
    % RAM use increases exponentially with ncon)
    perminds = nuniquereturns(@()logical(rand([1,ncon],'single')>.5),...
        nperms-1,[1,ncon]);
    % insert true perm (all false) on top
    perminds = vertcat(false([1,ncon]),perminds{:});
end

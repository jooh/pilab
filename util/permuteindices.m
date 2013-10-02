% generate nperms permutations of the indices in 1:ncon. The returned
% matrix perminds has the 'true' permutation (1:ncon) in the first row.
% perminds = permuteindices(ncon,nperms)
function perminds = permuteindices(ncon,nperms)

npossible = factorial(ncon);
if nperms > npossible
    warning(...
    'requested %d permutations but only %d are possible',nperms,npossible);
    nperms = npossible;
end

% find row/column indices for each permutation
if ncon<=10
    % extra operation to put the unpermuted entry first
    perminds = (ncon+1) - perms(1:ncon);
    % strip out original (to be reinserted later)
    perminds(1,:) = [];
    % randomise order
    perminds = perminds(randperm(npossible-1),:);
    % and restrict to nperms (keeping original as 1)
    perminds = perminds(1:nperms-1,:);
else
    % not enough memory to generate all perms so just draw randoms
    perminds = nuniquereturns(@()randperm(ncon),nperms-1,[1 ncon]);
    perminds = vertcat(perminds{:});
end

% return original perm
perminds = [1:ncon; perminds];

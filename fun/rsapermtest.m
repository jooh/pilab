% Permute rows and columns of rdm a and use the resulting null distribution
% to compute a p value for its rank relationship with each of the rdms in
% b.  Uses parfor.
%
% [r,p,nulldist] = rsapermtest(a,b,[nperms])
%
% INPUTS
% a: a single predictor RDM (struct / mat / vector form)
% b: one or many data RDMs
% nperms: defaults to 1000
%
% OUTPUTS
% r: 1 by n matrix of correlations between a and each entry in b
% p: permutation p values (one-tailed)
% nulldist: nperms by n matrix of r for each permutation (true is first)
function [r,p,nulldist] = rsapermtest(a,b,nperms)

if ieNotDefined('nperms')
    nperms = 1000;
end

% analyse how many true conditions exist (ie, rows where not all the
% off-diagonal entries are NaN). It doesn't make sense to permute empty
% rows/columns.
amat = asrdmmat(a);
% set diagonal to NaN to make the next test easier
amat(logical(eye(size(amat,1)))) = NaN;
ok = ~all(isnan(amat),1);
if ~all(ok)
    % shrink a and b down to size
    a = asrdmmat(a);
    a = a(ok,ok,:);
    b = asrdmmat(b);
    b = b(ok,ok,:);
end

% rank transform for speed (no need to recompute on each iteration)
a = tiedrank(asrdmvec(a));
b = tiedrank(asrdmvec(b));
ncon = npairs2n(length(a));
ndata = size(b,2);

npossible = factorial(ncon);
assert(nperms<=npossible,...
    'requested %d permutations but only %d are possible',nperms,npossible);

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
    done = 0;
    while ~done
        perminds = cell2mat(arrayfun(@randperm,ones(nperms-1,1)*ncon,...
            'uniformoutput',false));
        % it is usually exceedingly unlikely, but just out of paranoia we
        % will make sure that we don't draw the same perm twice, and that
        % the perms don't include the original
        if (size(unique(perminds,'rows'),1)==(nperms-1)) && ...
                isempty(intersect(perminds,1:ncon,'rows'))
            done = 1;
        end
    end
end
% return original perm
perminds = [1:ncon; perminds];

% preallocate
nulldist = NaN([nperms ndata]);

% compute true statistic
% (nb we use pearsonvec on ranks which is equivalent to spearman but faster
% since rank transform is only done once for b)
r = pearsonvec(a,b);
% return a to mat form for permuting
a = vec2rdm(a);
nulldist(1,:) = r;
% predefine column indices to prevent parfor problems
inds = 1:ndata;
% todo: parfor
parfor p = 2:nperms
    % vector form of permuted a
    aperm = rdm2vec(a(perminds(p,:),perminds(p,:)));
    nulldist(p,inds) = pearsonvec(aperm,b);
end

% get p for each RDM
p = sum(nulldist >= repmat(r,[nperms 1]),1) ./ nperms;

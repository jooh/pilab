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
% test for asymmetries in nans
assert(all(ok==~all(isnan(amat),2)'),['detected different nans in ' ... 
    'rows and columns of a - asymmetrical rdm?']);
% enforce matrix form
a = asrdmmat(a);
b = asrdmmat(b);
% remove nan'ed rows/columns
if ~all(ok)
    % shrink a and b down to size
    assert(size(a,1)==size(b,1),'different size rdms detected');
    a = a(ok,ok,:);
    b = b(ok,ok,:);
end
ncon = size(a,1);

% by default we use the stock rsa conversions between square and vector
% forms (basically squareform)
makevector = @(rdm) rdm2vec(rdm);
makerdmmat = @(rdvec) vec2rdm(rdvec);
splitdatamode = false;

% detect remaining NaNs in model
if any(isnan(a(:)))
    if issplitdatardm(a)
        splitdatamode = true;
        fprintf(['split data RDM detected, switching to row and column '...
            'permutation test\n']);
        % split data RDM test mode
        % reduce model to split RDM part
        % update ncon
        ncon = size(a,1)/2;
        % number of pairs for split data case (not nchoosek(x,2))
        npairs = ncon^2;
        rowinds = ncon+1:ncon*2;
        colinds = 1:ncon;
        a = a(rowinds,colinds);
        % also data
        b = b(rowinds,colinds,:);
        % switch vectorisation function to use reshape (to include
        % diagonal)
        makevector = @(rdm) reshape(rdm,[npairs size(rdm,3)]);
        makerdmmat = @(rdvec) reshape(rdvec,[ncon ncon size(rdvec,2)]);
        % switch permutation mode to row/columns
    else
        % some other NaNs in model - probably user error
        error('non-row/column NaNs in non-split data model RDM');
    end
end

% both back to vector - either unique off-diagonals or everything if split
% data
b = makevector(b);
a = makevector(a);
assert(~any(isnan(b(:))),'nans detected in data');
% rank transform for speed (no need to recompute on each iteration)
a = tiedrank(a);
b = tiedrank(b);
assert(isequal(size(a,1),size(b,1)),'different size rdms detected');
nfeatures = size(b,2);

npossible = factorial(ncon);
if npossible < nperms
    warning(['requested %d permutations but only %d are possible - ' ...
        'skipping permutation test.'],nperms,npossible);
    nperms = 1;
end

permrows = permuteindices(ncon,nperms);
if splitdatamode
    % for split data we permute rows and columns separately
    done = false;
    niter = 0;
    while ~done
        permcols = permuteindices(ncon,nperms);
        % toss symmetrical row/column permutations 
        if nperms==1 || ~any(all(permcols(2:end,:)==permrows(2:end,:),2))
            done = true;
        end
        niter = niter + 1;
        assert(niter < 1e5,'no valid row/column permutations found');
    end
else
    % for non-split data we permute rows and columns in concert to preserve
    % diagonal symmetry
    permcols = permrows;
end

% preallocate
nulldist = NaN([nperms nfeatures]);

% compute true statistic
% (nb we use pearsonvec on ranks which is equivalent to spearman but faster
% since rank transform is only done once for b)
r = pearsonvec(a,b);
nulldist(1,:) = r;
% return a to mat form for permuting
a = makerdmmat(a);
% predefine column indices to prevent parfor problems
inds = 1:nfeatures;
parfor p = 2:nperms
    % vector form of permuted a
    aperm = makevector(a(permrows(p,:),permcols(p,:)));
    nulldist(p,inds) = pearsonvec(aperm,b);
end

% get p for each RDM
if nperms>1
    p = sum(nulldist >= repmat(r,[nperms 1]),1) ./ nperms;
else
    p = NaN(size(r));
    nulldist(:) = NaN;
end

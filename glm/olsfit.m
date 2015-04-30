% Fit a design matrix X to responses in Y using ordinary least squares.
% Crash if the design matrix is rank deficient.
%
% NB, stores X,Y and betas in an internal cache to speed up re-fits. This
% has memory implications.
%
% betas = olsfit(X,Y)
function betas = olsfit(X,Y)

persistent olsfit_cache_X
persistent olsfit_cache_Y
persistent olsfit_cache_betas

if ~iscell(Y)
    Y = {Y};
end

if isequal(olsfit_cache_X,X) && isequal(olsfit_cache_Y,Y)
    betas = olsfit_cache_betas;
    return;
end
lastwarn('')
% nb slightly different formulation from the stock X \ Y - this is faster
% for large Ys.
% also, mtimescell - save memory by fitting each entry in the cell array Y
% separately
betas = mtimescell((X'*X)\X',Y);
assert(isempty(lastwarn),'rank deficient fit');
olsfit_cache_X = X;
olsfit_cache_Y = Y;
olsfit_cache_betas = betas;

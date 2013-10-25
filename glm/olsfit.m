% Fit a design matrix X to responses in Y using ordinary least squares.
% Crash if the design matrix is rank deficient.
%
% NB, stores X,Y and betas in an internal cache to speed up re-fits. This
% has memory implications.
%
% betas = olsfit(X,Y)
function betas = olsfit(X,Y)

global olsfit_cache_X
global olsfit_cache_Y
global olsfit_cache_betas

if isequal(olsfit_cache_X,X) && isequal(olsfit_cache_Y,Y)
    betas = olsfit_cache_betas;
    return;
end
lastwarn('')
betas = X \ Y;
assert(isempty(lastwarn),'rank deficient fit');
olsfit_cache_X = X;
olsfit_cache_Y = Y;
olsfit_cache_betas = betas;

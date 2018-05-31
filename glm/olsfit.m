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

% this approach is from Kendrick's mtimescell
Xt = (X'*X)\X';
betas = zeros([size(X,2),size(Y{1},2)],'like',Y{1});
cnt = 0;
for q=1:length(Y)
    betas = betas + Xt(:,cnt + (1:size(Y{q},1))) * Y{q};
    cnt = cnt + size(Y{q},1);
end
assert(isempty(lastwarn),'rank deficient fit');
olsfit_cache_X = X;
olsfit_cache_Y = Y;
olsfit_cache_betas = betas;

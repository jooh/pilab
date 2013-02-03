% Fit a design matrix X to responses in Y using ordinary least squares.
% Crash if the design matrix is rank deficient.
% betas = olsfit(X,Y)
function betas = olsfit(X,Y)

lastwarn('')
betas = X \ Y;
assert(isempty(lastwarn),'rank deficient fit');

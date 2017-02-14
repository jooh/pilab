% one-liner function to calculate the variance inflation index (VIF) for the
% design matrix X (samples by regressors).
%
% VIF is a ratio that summarises the extent to which each predictor can be
% explained as a linear combination of the remaining predictors. VIF=1 for
% orthogonal designs, VIF>5 is considered problematic in fMRI designs (Mumford),
% VIF=inf for rank deficient designs, but you in practice will likely obtain huge
% values instead from this function due to the instability of Matlab inv.
%
% v = vif(X)
function v = vif(X)

% convenient vectorised VIF calculation
v = diag(inv(corrcoef(X)));

% an equivalent but more didactic algorithm would be
% nreg = size(X,2);
% for n = 1:nreg
    %b = X(:,setdiff(1:nreg,n)) \ X(:,n);
    %yh = X(:,setdiff(1:nreg,n)) * b;
    %vifalt(n) = 1 / (1-corr(X(:,n),yh).^2);
%end

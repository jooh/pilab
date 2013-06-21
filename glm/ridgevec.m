% Vectorised ridge regression. Almost entirely snatched from Matlab's
% built-in ridge function but without the arbitrary restrction to only fit
% a single y per call. As a tradeoff, you can only fit one k per call. For
% general documentation, see ridge.
%
% y: data (nsamples by nfeatures)
% X: design matrix (nsamples by nregressors).
% k: ridge parameter (multiplied by n internally). 0 for OLS fit.
% unscale: (default 0) rescale fit with Z-scored design matrix to original
%   units (with constant term - gets inserted in last column)
%
% Changes to Matlab's default ridge behaviour:
% a) if unscaling is off, no constant term is added
% b) the ridge parameter (k) is scaled by n (Draper & Smith, 1998). This
% tends to mean that you don't have to go to extremely large k to see a
% divergence from the k=0 case.
%
% b = ridgevec(y,X,k,[unscale])
function b = ridgevec(y,X,k,unscale)

if ieNotDefined('unscale')
    unscale = 0;
end

[n,p] = size(X);
% Draper & Smith scaling
k = k*n;

[n1,collhs] = size(y);
if n~=n1, 
    error('ridgevec:InputSizeMismatch',...
          'The number of rows in Y must equal the number of rows in X.'); 
end 

% Remove any missing values
wasnan = ((any(isnan(y),2)) | any(isnan(X),2));
if (any(wasnan))
   y(wasnan,:) = [];
   X(wasnan,:) = [];
end

% Normalize the columns of X to mean zero, and standard deviation one.
mx = mean(X,1);
stdx = std(X,0,1);
% catch bad user behaviour
assert(~any(stdx==0),...
    'Constant detected. Fit without it and use unscale==1 instead.');
idx = find(abs(stdx) < sqrt(eps(class(stdx)))); 
if any(idx)
  stdx(idx) = 1;
end

MX = mx(ones(n,1),:);
STDX = stdx(ones(n,1),:);
Z = (X - MX) ./ STDX;
if any(idx)
  Z(:,idx) = 1;
end

% Compute the ridge coefficient estimates using the technique of
% adding pseudo observations having y=0 and X'X = k*I.
pseudo = sqrt(k) * eye(p);
Zplus  = [Z;pseudo];
yplus  = [y;zeros(p,collhs)];

% and now the b is...
b = olsfit(Zplus,yplus);

% Return fit with Z-scored DM to fit from unscaled DM with constant term?
if unscale
    % each beta divided by std of design matrix (equivalent to fitting with
    % constant that is not penalised)
    b = b ./ repmat(stdx',[1 collhs]);
    % mean of Y - mean of X (equivalent to constant fit in design matrix)
end

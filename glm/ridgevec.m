% Vectorised ridge regression. Almost entirely snatched from Matlab's
% built-in ridge function but without the arbitrary restrction to only fit
% a single y per call. As a tradeoff, you can only fit one k per call. For
% documentation, see ridge.
%
% Note that this function sometimes produces results that are near but not
% absolutely identical to the Matlab built-in. This problem originates in
% strange behaviour in Matlab's backslash operator, where the first
% parameter estimate from  dm\y does not always equal the only estimate
% from dm\y(:,1). Anyway, any deviations should be extremely subtle.
%
% b = ridgevec(y,X,k)
function b = ridgevec(y,X,k)

% from builtin ridge
[n,p] = size(X);

[n1,collhs] = size(y);
if n~=n1, 
    error('stats:ridge:InputSizeMismatch',...
          'The number of rows in Y must equal the number of rows in X.'); 
end 

% Remove any missing values
wasnan = ((any(isnan(y),2)) | any(isnan(X),2));
if (any(wasnan))
   y(wasnan,:) = [];
   X(wasnan,:) = [];
   n = length(y);
end

% Normalize the columns of X to mean zero, and standard deviation one.
mx = mean(X);
stdx = std(X,0,1);
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
b = Zplus\yplus;

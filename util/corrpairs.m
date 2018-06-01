% return a vector of the paired correlations between the columns in a and
% b (over rows).
%
% r = corrpairs(a,b)
function r = corrpairs(a,b)

r = dot(unitlen(bsxfun(@minus,a,mean(a))),unitlen(bsxfun(@minus,b,...
    mean(b))));

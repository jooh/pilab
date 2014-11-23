% return a vector of the paired correlations between the columns in a and
% b (over rows).
%
% r = corrpairs(a,b)
function r = corrpairs(a,b)

r = dot(unitlen(zeromean(a,1)),unitlen(zeromean(b,1)));

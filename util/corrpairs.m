% return a vector of the paired correlations between the columns in a and
% b (over rows).
%
% r = corrpairs(a,b)
function r = corrpairs(a,b)

r = dot(unitlength(zeromean(a,1),1),unitlength(zeromean(b,1),1));

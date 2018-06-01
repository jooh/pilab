% project out (ie subtract the fitted contribution of) covariate x on data
% y. Note that this is faster than premultiplying a projection matrix (ie, 
% y = projectionmatrix(c) * y) when working with small matrices, but the
% difference then reverses quite dramatically. Worth testing with your own
% code which is faster.
%
% y = projectout(y,c)
function y = projectout(y,c)

% data - fitted response
y = y - (c * (c \ y));

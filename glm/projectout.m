% project out (ie subtract the fitted contribution of) covariate x on data
% y. 
%
% y = projectout(y,c)
function y = projectout(y,c)

% data - fitted response
y = y - (c * (c \ y));

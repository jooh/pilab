% Get R^2 for a set of fitted responses. (Yfit = X * betas).
% R = rsquare(Yfit,Y)
function R = rsquare(Yfit,Y)

SSerr = sum((Yfit-Y).^2,1);
SStot = (size(Y,1)-1) .* var(Y);
R = 1 - (SSerr ./ SStot);

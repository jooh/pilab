% Get mean square error for a set of fitted responses. (Yfit = X * betas).
% mse = mserr(Yfit,Y)
function mse = mserr(Yfit,Y)

mse = mean((Yfit-Y).^2);

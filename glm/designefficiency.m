% return the design efficiency for each contrast (row) in c given the
% (convolved, filtered) design matrix X. 
%
% Given a fixed noise level, design efficiency is monotonically (but not
% linearly) related to the standard error of the contrast estimate. So design
% matrices with higher efficiency for a given contrast are expected to produce
% more stable estimates.
%
% See also: vif
%
% 20171003 J Carlin
%
% eff = designefficiency(X,c)
function eff = designefficiency(X,c)

assert(size(X,2) == size(c,2),...
    'contrast vector c must have same number of columns as the design matrix X')

eff = diag(c * inv(X'*X) * c');

% Kendall's tau-a. Adapted from original by Niko.
%
% taua = kendall_a(a,b)
function taua= kendall_a(a,b)

n=size(a,1);

%% compute Kendall rank correlation coefficient tau-a
K = 0;
for k = 1:n-1
    pairRelations_a=sign(a(k)-a(k+1:n));
    pairRelations_b=sign(b(k)-b(k+1:n));
    K = K + sum(pairRelations_a.*pairRelations_b);
end
taua=K/(n*(n-1)/2); % normalise by the total number of pairs 

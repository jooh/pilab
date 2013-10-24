% return true if the entered matrix is a valid rdm (isrdm test), and has at
% most 1 dissimilarity in diag(rdm,-1).
%
% yes = issplitdatardm(rdm)
function yes = issplitdatardm(rdm)

yes = true;

if ~isrdm(rdm)
    yes = false;
    return;
end

nrdm = size(rdm,3);
for n = 1:nrdm
    if sum(~isnan(diag(rdm(:,:,n),-1)))~=1
        yes = false;
        return;
    end
end

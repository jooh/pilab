% return true if the entered matrix is square, has 0 diagonals and is
% symmetric.
% yes = isrdm(rdm)
function yes = isrdm(rdm)

yes = true;

[r,c,z] = size(rdm);
% square matrix at the very least
if r ~= c
    yes = false;
    return
end

% 0 diagonal for all RDMs (stacked in 3rd dim)
diagonals = rdm(repmat(logical(eye(r)),[1 1 z]));
if ~all(diagonals==0) && ~all(diagonals==1)
    yes = false;
    return
end

% symmetry
symtest = rdm == permute(rdm,[2 1 3]);
% need to exclude NaNs because oddly, NaN ~= NaN in Matlab
nans = isnan(rdm);
if ~all(symtest(~isnan(rdm)))
    yes = false;
end

% return true if the entered matrix is square, has 0 or 1 diagonals and is
% symmetric. Also supports RDM structs (must have fields name and RDM, and
% the RDM field must also pass the isrdm test - we recurse so nested RDM
% groups should work).
%
% yes = isrdm(rdm)
function yes = isrdm(rdm)

yes = true;

if isstruct(rdm)
    if ~isfield(rdm,'RDM')
        yes = false;
        return;
    end
    if ~isfield(rdm,'name')
        yes = false;
        return;
    end
    if ~isrdm(rdm(1).RDM)
        yes = false;
        return;
    end
else
    [r,c,z] = size(rdm);
    % square matrix at the very least
    if r ~= c
        yes = false;
        return
    end

    if r<2
        % too small
        yes = false;
        return;
    end

    % 0 diagonal for all RDMs (stacked in 3rd dim)
    diagonals = rdm(diagind([r c z]));
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
end

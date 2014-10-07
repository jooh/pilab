% convenience function for setting the diagonal entries in the input RDMs
% (matrix form, possibly stacked in 3rd dim) to 0.
%
% rdm = zerodiagonal(rdm)
function rdm = zerodiagonal(rdm)

rdm(diagind(size(rdm))) = 0;

% Convert a set of vector or struct RDMs to 3D matrix form. Do nothing if
% the RDMs are already in this format. Note that matrix RDMs MUST have
% zeros on the diagonal to prevent unexpected results. Use isrdm to test
% this.
% rdmmat = asrdmmat(in)
function rdmmat = asrdmmat(rdms)

if isnumeric(rdms)
    if isrdm(rdms)
        rdmmat = rdms;
    else
        % vectorised rdms then?
        assert(size(rdms,3)==1,'vectorised RDMs must be 2D')
        rdmmat = vec2rdm(rdms);
    end
elseif isstruct(rdms)
    % assume RSA style struct array
    n = size(rdms(1).RDM,1);
    rdmmat = reshape([rdms.RDM],[n n length(rdms)]);
else
    error('received rdms in unknown form')
end

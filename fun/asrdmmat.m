% Convert a set of vector or struct RDMs to 3D matrix form. Do nothing if
% the RDMs are already in this format
% rdmmat = asrdmmat(in)
function rdmmat = asrdmmat(rdms)

if isnumeric(rdms)
    [r,c,z] = size(rdms);
    if r==c && all(diag(rdms(:,:,1))==0)
        % assume already rdmmat
        rdmmat = rdms;
    else
        % vectorised rdms then?
        assert(z==1,'vectorised RDMs must be 2D')
        rdmmat = vec2rdm(rdms);
    end
elseif isstruct(rdms)
    % assume RSA style struct array
    n = size(rdms(1).RDM,1);
    rdmmat = reshape([rdms.RDM],[n n length(rdms)]);
else
    error('received rdms in unknown form')
end

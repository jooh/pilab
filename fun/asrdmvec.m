% Convert a set of matrix or struct RDMs to vector form. Do nothing if
% the RDMs are already in this format
% rdmvec = asrdmvec(in)
function rdmvec = asrdmvec(rdms)

if isnumeric(rdms)
    [r,c,z] = size(rdms);
    if r==c && all(diag(rdms(:,:,1))==0)
        % assume rdmmat
        rdmvec = rdm2vec(rdms);
    else
        % vectorised rdms then?
        assert(z==1,'vectorised RDMs must be 2D')
        rdmvec = rdms;
    end
elseif isstruct(rdms)
    % assume RSA style struct array
    n = size(rdms(1).RDM,1);
    rdmvec = rdm2vec(reshape([rdms.RDM],[n n length(rdms)]));
else
    error('received rdms in unknown form')
end

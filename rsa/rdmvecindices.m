% return logical indices corresponding to particular rows and columns in
% a symmetrical, square RDM matrix. conind is either logical or a
% row/column and can contain one or multiple entries.
% inds = rdmvecindices(conind,ndis)
function inds = rdmvecindices(conind,ndis)

if islogical(conind)
    conind = find(conind);
end
ncon = npairs2n(ndis);
assert(ncon==round(ncon),'ndis: %d invalid size for a vectorised RDM',ndis);
nind = length(conind);
assert(all(conind > 0 & conind <= ncon),...
    'conind: %d is 0 or too large for %d x %d RDM',...
    conind,ncon,ncon);

% this is a bit hacky but very safe
d = false([ncon ncon]);
d(conind,:) = true;
d(:,conind) = true;
d(logical(eye(ncon))) = false;
inds = squareform(d);

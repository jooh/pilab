% firdm = rdm2fir(rdm)
function firdm = rdm2fir(rdm)

rdvec = asrdmvec(rdm);
[ndis,nrdm] = size(rdvec);
assert(nrdm==1,'only one input RDM is supported');

if isstruct(rdm)
    name = rdm.name;
else
    name = 'firdm';
end

u = unique(rdvec);
u(isnan(u)) = [];
u(u==0) = [];
nfir = numel(u);

for f = 1:nfir
    newvec = zeros([ndis 1],class(rdvec));
    newvec(rdvec==u(f)) = 1;
    firdm(f) = struct('name',sprintf('%s_f%02.0f',name,u(f)),...
        'RDM',vec2rdm(newvec));
end

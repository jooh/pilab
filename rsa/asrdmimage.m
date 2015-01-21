% convert some rdm to image format for visualisation. Returns one cell
% array per RDM if multiple rdms are entered. If the input is a struct with
% an im field we return this instead.
%
% im = asrdmimage(rd,[cmap])
function im = asrdmimage(rd,cmap)

if isstruct(rd) && isfield(rd,'im')
    im = {rd.im};
end

if ieNotDefined('cmap')
    cmap = [];
end

rdm = asrdmmat(rd);
for r = 1:size(rdm,3)
    im{r} = intensity2rgb(rdm(:,:,r),cmap);
    im{r}(diagind(size(im{r}))) = 1;
end

if numel(im)==1
    im = im{1};
end

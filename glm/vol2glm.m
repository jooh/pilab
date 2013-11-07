% Convert a pair of Volume-derived instances (design and data) to a
% GLM-derived object array instance by iterating over the unique chunks.
%
% inputs:
% designvol: data field contains design matrix.
% datavol: data field contains samples.
% glmclass: (default GLM) type of GLM class to initialise
% glmargs: (default []) other arguments to pass on to GLM class
%
% glm = vol2glm(designvol,datavol,glmclass,glmargs)
function glm = vol2glm(designvol,datavol,glmclass,varargin)

chunks = designvol.desc.samples.unique.chunks;
labels = designvol.meta.features.labels;
for c = 1:designvol.desc.samples.nunique.chunks
    thischunk = chunks(c);
    glm(c) = feval(glmclass,designvol.data(...
        designvol.meta.samples.chunks==thischunk,:),datavol.data(...
        datavol.meta.samples.chunks==thischunk,:),varargin{:});
end

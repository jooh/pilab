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

if ieNotDefined('glmclass')
    glmclass = 'GLM';
end

chunks = designvol.desc.samples.unique.chunks;
glm = arrayfun(@(c) feval(glmclass,...
    designvol.data(designvol.meta.samples.chunks==c,:),...
    datavol.data(datavol.meta.samples.chunks==c,:),varargin{:}),chunks,...
    'uniformoutput',false);
glm = vertcat(glm{:});

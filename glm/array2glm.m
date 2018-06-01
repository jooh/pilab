% Convert a set of arrays to a GLM-derived object array instance by
% iterating over unique entries in chunks. This is potentially faster than
% vol2glm.
%
% inputs:
% designmat: nsamples by npredictors
% datamat: nsamples by nfeatures
% chunks: nsamples by 1
% glmclass: (default GLM) type of GLM class to initialise
% glmargs: (default []) other arguments to pass on to GLM class
%
% glm = array2glm(designmat,datamat,chunks,glmclass,[glmargs])
function glm = array2glm(designmat,datamat,chunks,glmclass,varargin)

uc = unique(chunks);
% it is tempting to attempt to pre-allocate the GLM instance here but this
% would make very strong assumptions on chunks (identical shape of all
% data). 
for c = 1:numel(uc)
    glm(c) = feval(glmclass,designmat(chunks==uc(c),:),...
        datamat(chunks==uc(c),:),varargin{:});
end

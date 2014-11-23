% Create a GLM-derived instance (by default RidgeRSA) from a set of model
% RDMs (any pilab RSA format or a Volume instance) and a cell array of
% disvols. See also vol2glm. We assume that you want the same dissimilarity
% design matrix for all splits.
%
% glm = disvol2glm(model,disvols,[glmclass],[varargin]);
function glm = disvol2glm(model,disvols,glmclass,varargin);

if ieNotDefined('glmclass')
    glmclass = 'RankRSA';
end

if isa(model,'Volume')
    model = model.data;
end
% otherwise assume you entered an appropriately structured set pilab RDMs.

glm = arrayfun(@(s) feval(glmclass,model,disvols{s}.data,varargin{:}),...
    1:length(disvols),'uniformoutput',false);
glm = vertcat(glm{:});

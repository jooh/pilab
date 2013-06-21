% convert a SPM struct to a GLM instance. Data are loaded in 2d masked
% form. Non-task regressors in SPM are projected out of the data and design
% matrix for each run. 
%
% SPM - SPM.mat struct
% mask - 3D matrix of same shape as EPI (e.g. the mask.hdr generated during
%   model estimation)
% glmclass - GLM-subclass flavour (default GLM)
% glmargs - any additional arguments are passed on to glmclass
%
% glm = spm2glm(SPM,mask,[glmclass],[glmargs]);
function glm = spm2glm(SPM,mask,glmclass,varargin);

if ieNotDefined('glmclass')
    glmclass = 'GLM';
end

% load data
fprintf('loading masked EPI volumes...\n')
[data,mask] = loadmaskedvolumes(SPM.xY.P,mask);

% configure design matrices
x = SPM.xX.X;
[nvol,nreg] = size(SPM.xX.X);
nrun = length(SPM.Sess);
for r = 1:nrun
    fprintf('configuring instance for run %d of %d...\n',r,nrun);
    % separate task regressors from nuisance covariates
    taskreg = SPM.Sess(r).col(1:length(SPM.Sess(r).U));
    nuisancereg = setdiff(SPM.Sess(r).col,taskreg);
    runx = x(SPM.Sess(r).row,taskreg);
    nuisancex = x(SPM.Sess(r).row,nuisancereg);
    if ~any(all(nuisancex==1,1))
        % re-insert run constant
        nuisancex = [nuisancex ones(size(runx,1),1)];
    end
    % project out nuisance covariates
    nmat = projectionmatrix(nuisancex);
    glm(r) = feval(glmclass,nmat*runx,nmat*data(SPM.Sess(r).row,:),...
        varargin{:});
end

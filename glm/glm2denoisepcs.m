function [noisepcs,result] = glm2denoisepcs(model,trainind,varargin)

assert(isa(model,'CovariateGLM'),['input glm instance must be ' ...
    'CovariateGLM or sub-class thereof']);

numpcstotry = [];
getArgs(varargin,{'brainthresh',[99 .5],'brainR2',0,'brainexclude',...
    false([1,model(1).nfeatures]),'numpcstotry',20});

% make sure we are using all the available covariates in each run
maxn = num2cell([model.ncovariates]);
[model.ncovtouse] = deal(maxn{:});

% identify noise pool based on training part of instance vector
result.r2before = cvpredictionrun(model(trainind),'predictY','rsquare');
result.meanepi = mean(vertcat(model(trainind).data));
result.thresh = prctile(result.meanepi(:),brainthresh(1)) * brainthresh(2);
result.bright = result.meanepi > result.thresh;
noisepool = result.bright & result.r2before<brainR2 & ~brainexclude;

% but generate noise PCs for all instances
datac = getdatac(model);
parfor n = 1:numel(model)
    % get noisepool time series
    thisdata = unitlen(datac{n}(:,noisepool));
    % obtain the desired number of PCs
    [noisepcs{n},~,~] = svds(double(thisdata*thisdata'),numpcstotry);
    % scale to st dev 1 and convert class back to input
    noisepcs{n} = cast(bsxfun(@rdivide,noisepcs{n},...
        std(noisepcs{n},[],1)),'like',thisdata);
end
result.noisepool = noisepool;

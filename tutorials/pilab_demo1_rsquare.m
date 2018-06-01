% demonstrate how to use pilab to fit and cross-validate a simple fMRI model.
% Raw matlab code version of the ipython notebook, which contains more
% documentation and explanation.
function pilab_demo1_rsquare()

% shuffle randomisation
rng('shuffle');
% control variables for the simulation
% number of volumes (TRs)
n = 200;
% frame period - one volume every tr s
tr = 2;
% number of onsets to include in each run
ntrial = 20;
ncon = 4;
nvox = 50;
nrun = 2;
% true parameters - same across runs. +1 to model variations in constant
b = rand(ncon+1,nvox);

designvol = cell(nrun,1);
datavol = cell(nrun,1);
for r = 1:nrun
    % different trial order in each run
    conind = ceil(rand(ntrial,1) * ncon);
    onsets = rand(ntrial,1) * n;
    X = [convolveonsets(onsets,conind,tr,n,'dur',1) ones(n,1)];
    % construct data as some mixture of gaussian noise and the true parameters
    % scaled by each predictor
    data = randn(n,nvox)*.3*max(X(:)) + X*b;
    % build Volume instances for convolved design and data
    metasamp = struct('chunks',ones(n,1)*r);
    designvol{r} = Volume(X,'frameperiod',tr,'metasamples',metasamp,...
        'metafeatures',struct('labels',{{'a','b','c','d','constant'}}));
    datavol{r} = Volume(data,'frameperiod',tr,'metasamples',metasamp);
end
% operator overloading means we can use standard concatenation syntax to combine
% the Volume instances over runs
datavol = cat(1,datavol{:});
designvol = cat(1,designvol{:});
datavol

%plot -s 700,700 -r 300 -f svg
ax = subplot(2,1,1);
imagesc(datavol.data');
title('data');
ylabel('voxels');
% hide the x axis. Don't you wish Matlab provided a more elegant way of doing this?
set(gca,'xcolor',[1 1 1]);
ax(2) = subplot(2,1,2);
% plot the data in second units by scaling frameperiod
ph = plot((1:designvol.nsamples)*designvol.frameperiod,designvol.data);
xlabel('time (s)');
title('model');
lh = legend(ph,designvol.meta.features.labels,'location','southoutside');
set([lh ax],'box','off');

[filtereddatavol,filtereddesignvol] = preprocessvols(datavol,designvol,...
    'covariatedeg','','ignorelabels','constant');
model = vol2glm(filtereddesignvol,filtereddatavol);
% notice that the model only contains 4 regressors as the constant has been projected out
model(1)

%plot -s 700,350 -r 300 -f svg
% cross-validated prediction - use model(2) fit and model(1) design matrix
% Notice that the GLM is an object array class, so can be indexed like a struct
% to run particular methods on say a train and a test split of the data.
prediction = predictY(model(2),model(1).X);
% evaluate on model(1) data
r2 = rsquare(model(1),prediction)
% visualise fit
ph = plot([prediction(:,1) model(1).data(:,1)],'.-');
lh = legend(ph,{'prediction','data'},'location','best');
set([gca,lh],'box','off');

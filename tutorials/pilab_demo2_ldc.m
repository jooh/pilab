% demonstrate how to use pilab to estimate linear discriminant contrast, a
% pseudo-distance metric. Raw matlab code version of the ipython notebook, which
% contains more documentation and explanation.
function pilab_demo2_ldc()

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
    X = convolveonsets(onsets,conind,tr,n,'dur',1);
    % scale design matrix to peak at 1
    X = X / max(X(:));
    % add a constant at the end
    X = [X ones(n,1)];
    % construct data as some mixture of gaussian noise and the true parameters
    % scaled by each predictor
    data = randn(n,nvox)*.3 + X*b;
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

[filtereddatavol,filtereddesignvol] = preprocessvols(datavol,designvol,...
    'covariatedeg','','ignorelabels','constant');
model = vol2glm(filtereddesignvol,filtereddatavol);

% cross-validated linear discriminant contrast (train on 1, test on 2)
m = infoc(model(2),discriminant(model(1),[1 -1 0 0]),[1 -1 0 0])

% OLS fit of data from run 1
est1 = model(1).X \ model(1).data;
% contrast between the conditions of interest
con1 = [1 -1 0 0] * est1;
% fitted time course
predict1 = model(1).X * est1;
% residual time courses
res1 = model(1).data - predict1;
% (sparse) covariance of the residuals - see Misaki et al., 2010, NeuroImage. You could use the
% built-in cov function with very similar results.
cov1 = covdiag(res1);
% now the linear discriminant weights vector is just the contrast multiplied by the inverse of the
% covariance estimator - note that so far, nothing we have done differs from bog-standard LDA.
w1 = con1 / cov1;
% additionally we scale the weights to unit length, which turns out to be useful for expressing
% the cross-validated discriminant in the units of the test data.
w1 = unitlen(w1);
% now we get to the part where the discriminant contrast differs from conventional classification,
% where you would attempt to predict the label for individual observations in the test data. Instead,
% our discriminant contrast metric can be thought of as the distance of the (average) test data from
% the hyperplane defined by w1.
% OLS fit of data from run 2
est2 = model(2).X \ model(2).data;
% contrast between conditions of interest in run 2
con2 = [1 -1 0 0] * est2;
% and now the linear discriminant contrast is just the dot product between the weights vector
% and the contrast estimate from the test data. Intuitively, consistent signed values will sum
% and inconsistent signs will subtract from the final statistic.
m = con2 * w1'

% let's first build an additional region of interest Volume, which stores the valid
% features for each ROI
roivol = Volume([ones(1,nvox); ones(1,floor(nvox/2)) zeros(1,ceil(nvox/2)); ...
    zeros(1,floor(nvox/2)) ones(1,ceil(nvox/2))],'metasamples',struct(...
    'names',{{'allfeatures','firsthalf','secondhalf'}}));
% now we can just call this high-level analysis function to obtain the
% cross-validated linear discriminant contrast for each pair of conditions, and for
% each region
disvol = roidata2rdmvol_lindisc_batch(roivol,filtereddesignvol,filtereddatavol);

%plot -s 700,700 -r 300 -f svg
fh = gcf;
% round to 1 decimal place
maxd = reduceprecision(max(abs(disvol.data(:))),1);
[fh,ax,intmap,cmap] = rdmfig(disvol.data,disvol.meta.features.names,[],fh,...
    'limits',[-maxd maxd]);
cbax = subplot(2,2,4);
c = colorbarbetter(cbax,intmap{end},cmap{end},'label','linear discriminant contrast','scale',.5);

% to get from infoc to infot, just call a different GLM class method
t = infot(model(2),discriminant(model(1),[1 -1 0 0]),[1 -1 0 0])

% we go from one feature per voxel to one feature per discriminant
disctime = model(2).data * w1';

%plot -s 700,700 -r 300 -f svg
subplot(3,1,1);
plot(disctime);
title('data');
subplot(3,1,2);
plot(model(2).X(:,1:2));
title('regressors of interest');
subplot(3,1,3);
plot(model(2).X(:,1:2) * [1 -1]');
title('contrast time course');

% parameter estimates for the test design matrix fitted to the discriminant time course
discest = model(2).X \ disctime;
% contrast estimate
disccon = [1 -1 0 0] * discest
